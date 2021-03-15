using DataStructures
using FixedPointNumbers

struct T2{T}
  " creator "
  i :: T

  " annihilator "
  j :: T
end

struct T4{T}
  " creator "
  i :: T
  j :: T

  " annihilator "
  k :: T
  l :: T
end

const Kink{T} = Union{T2{T},T4{T}}

" outer constructor method to construct a T2 kink, inferring the type parameter from the arguments "
Kink(i,j) = T2(i,j)
" outer constructor method to construct a T4 kink, inferring the type parameter from the arguments "
Kink(i,j,k,l) = T4(i,j,k,l)
""" outer constructor method to extract a kink from a pair where the second element is a kink.
    This is useful for automatic conversion when looping over SortedDict{S,Kink{T}} """
Kink(p::Pair{S,T} where {T<:Kink} where {S}) = p[2]# first substitute S, then T

const ImgTime = Fixed{Int64,60}

" multi-particle trajectory using single particle states with type T "
mutable struct Configuration{T}
  " set of orbitals occupied at τ=0 "
  occupations :: Set{T}
  " excitations, using τ as an index "
  kinks :: SortedDict{ImgTime, Kink{T}}
end

" outer constructor method for a configuration with occupations given by o and kinks given by k. k can be anything from which a SortedDict can be constructed from. "
Configuration(o::Set{T}, k) where {T} = Configuration(o, SortedDict{ImgTime,Kink{T}}(Base.Forward, k) )
" outer constructor method for a configuration with occupations given by o and no kinks. "
Configuration(o::Set{T}) where {T} = Configuration(o, SortedDict{ImgTime,Kink{T}}(Base.Forward))
" outer constructor method for a configuriation with no occupations and kinks given by k. k can be anything from which a SortedDict can be constructed from. "
Configuration(k::SortedDict{ImgTime,<:Kink{T}}) where {T} = Configuration(Set{T}(), k)
" outer constructor method for a configuration with no occupations and kinks as given by varargs p... of Pair{ImgTime,<:Kink}, which are passed to the SortedDict constructor "
Configuration(p::Pair{ImgTime,<:Kink{T}}...) where {T} = Configuration(SortedDict{ImgTime,Kink{T}}(Base.Forward, p...))
" outer constructor method for a configuration with occupations given by o and kinks as given by varargs p... of Pair{ImgTime,<:Kink}, which are passed to the SortedDict constructor "
Configuration(o::Set{T}, p::Pair{ImgTime,<:Kink{T}}...) where {T} = Configuration(o, SortedDict{ImgTime,Kink{T}}(Base.Forward, p...))

" outer constructor method for empty Configurations{T} "
Configuration{T}() where T = Configuration(Set{T}())

" get single particle basis type "
basis(c::Configuration{T}) where T = T


abstract type Orbital end

#TODO Change the name of function to something more self explenatory for example "apply_T4"
" apply a T4 kink to a set of basis states "
function kink(o::Set{T}, κ::T4{T}) where T
  @assert ( in(κ.k, o) & in(κ.l, o) ) "Kink ($(κ.i),$(κ.j),$(κ.k),$(κ.l)) cannot be applied: one or two of the annihilators $(κ.k), $(κ.l) is not occupied. (Pauli-Principle)"
  @assert ( !in(κ.i, o) & !in(κ.j, o) ) "Kink ($(κ.i),$(κ.j),$(κ.k),$(κ.l)) cannot be applied: one or two of the creators $(κ.i), $(κ.j) is already occupied. (Pauli-Principle)"
  union(setdiff(o, Set([κ.k,κ.l])), Set([κ.i, κ.j]))
end

""" Apply a T4 kink to a set of basis states for a pair of a time and a T4-kink.
    This is useful for iteration of a SortedDict{ImgTime, T4{T}}."""
kink(o::Set{T}, κ::Pair{ImgTime,T4{T}}) where T = kink(o, κ[2])

" Apply a T4 kink in-place to a set of basis states. "
function kink!(o::Set{T}, κ::T4{T}) where T
  @assert (in(κ.k, o) & in(κ.l, o)) "Kink ($(κ.i),$(κ.j),$(κ.k),$(κ.l)) cannot be applied: one or two of the annihilators $(κ.k), $(κ.l) is not occupied. (Pauli-Principle)"
  @assert (!in(κ.i, o) & !in(κ.j, o)) "Kink ($(κ.i),$(κ.j),$(κ.k),$(κ.l)) cannot be applied: one or two of the creators $(κ.i), $(κ.j) is already occupied. (Pauli-Principle)"
  delete!(o, κ.k)
  delete!(o, κ.l)
  push!(o, κ.i)
  push!(o, κ.j)
end

" Return the occupied orbitals after applying all kinks to initial occupation. "
function occupations(o::Set{T}, kinks::SortedDict{ImgTime,Kink{T}}) :: Set{T} where {T}
  reduce(kink, kinks; init=o)
end

" Return the occupied orbitals to the right of τ ."
function occupations(c::Configuration, τ::ImgTime)
  occupations(c.occupations, filter(x -> x[1] <= τ, c.kinks))
end

" Return a list of all kinks that affect the given orbital. "
function get_kinks_of_orb(c::Configuration, orbital::Orbital)
  kinks_of_orb = SortedDict{ImgTime, Kink{<:Orbital}}()
  for (τ_kink,kink) in c.kinks
    if (kink.i == orbital) | (kink.j == orbital) | (kink.k == orbital) | (kink.l == orbital)
      kinks_of_orb[τ_kink] = kink
    end
  end
  return kinks_of_orb
end

""" Return a tuple of the imaginary times between 0 and 1 of the closest kink before τ and the closest kink after τ,
    which affect the given orbital. Has to be multiplied with β to get real imaginary times."""
function get_nearest_τ_affecting_orb(Configuration::Configuration, orbital::Orbital,τ::ImgTime)
  current_τ = 0
  kinks_of_orb = get_kinks_of_orb(Configuration, orbital)
  if length(kinks_of_orb) == 0
      return(nothing,nothing)
  end
  # TODO: binary search
  for (τ_kink,kink) in kinks_of_orb
      if τ_kink > τ
        if current_τ == 0
          return (first(last(kinks_of_orb)),τ_kink)
        else
          return (current_τ,τ_kink)
        end
      else
          # es soll nicht τ als Grenze zurückgegeben werden
          if τ_kink != τ
              current_τ = τ_kink
          end
      end
  end
  return (current_τ, first(first(kinks_of_orb)))
end



function get_τ_borders(Configuration::Configuration, orbitals::Set{<:Orbital},τ::ImgTime)
  if length(Configuration.kinks) == 0
      return(ImgTime(0),ImgTime(1))
  end
  #Initially we set τ right and τ left to the nearest Kinks left and right of τ.
  τ_left_semi_token  = searchsortedafter(Configuration.kinks, τ)
  τ_right_semi_token = searchsortedlast(Configuration.kinks, τ)
  if τ_left_semi_token == pastendsemitoken(Configuration.kinks)
      τ_left = first(first(Configuration.kinks))
  else
      τ_left = first(deref((Configuration.kinks, τ_left_semi_token)))
  end
  if τ_right_semi_token == beforestartsemitoken(Configuration.kinks)
      τ_right = first(last(Configuration.kinks))
  else
      τ_right = first(deref((Configuration.kinks, τ_right_semi_token)))
      if τ_right == τ
          τ_right_semi_token = regress((Configuration.kinks,τ_right_semi_token))
          if τ_right_semi_token == beforestartsemitoken(Configuration.kinks)
              τ_right = first(last(Configuration.kinks))
          else
              τ_right = first(deref((Configuration.kinks, τ_right_semi_token)))
          end
      end
  end
  # now search for the nearst τ's that do actually affect one of the orbitals
  non_interacting_orb_counter = 0
  for orb in orbitals
    @assert(τ_right != ImgTime(1))
    tupel = get_nearest_τ_affecting_orb(Configuration, orb, τ)
    if tupel[1] == nothing
        non_interacting_orb_counter += 1
    else
        # here we always have to check wether the given intervall extends over 1
        if τ_left < τ < tupel[1]
            nothing
        elseif tupel[1] < τ < τ_left
            τ_left = tupel[1]
        elseif τ_left < tupel[1]
          τ_left = tupel[1]
        end

        if tupel[2] < τ < τ_right
            nothing
        elseif τ_right < τ < tupel[2]
          τ_right = tupel[2]
        elseif tupel[2] < τ_right
            τ_right = tupel[2]
        end
    end
  end
  if non_interacting_orb_counter == length(orbitals)
      return (ImgTime(0),ImgTime(1))
  else
      return (τ_left,τ_right)
  end
end

" Return if an orbital is not affected by any kink. "
function is_non_interacting(Configuration::Configuration{T}, orbital::T) :: Bool where {T<:Orbital}
  for (τ_kink,kink) in Configuration.kinks
    if (kink.i == orbital) | (kink.j == orbital) | (kink.k == orbital) | (kink.l == orbital)
      return false
    end
  end
  return true
end


" Return if an orbital is not affected by any kink in the open interval (τ_first,τ_last). "
function is_non_interacting_in_interval(Configuration::Configuration{T}, orbital::T, τ_first::ImgTime, τ_last::ImgTime) :: Bool where {T<:Orbital}
  @assert τ_first != τ_last
  if τ_first < τ_last
      for (τ_kink,kink) in Configuration.kinks
        if (τ_kink <= τ_first) | (τ_kink >= τ_last)
            nothing
        elseif (kink.i == orbital) | (kink.j == orbital) | (kink.k == orbital) | (kink.l == orbital)
              return false
        end
      end
  else
      for (τ_kink,kink) in Configuration.kinks
        if ((τ_kink <= τ_first) & (τ_kink >= τ_last))
            nothing
        elseif (kink.i == orbital) | (kink.j == orbital) | (kink.k == orbital) | (kink.l == orbital)
              return false
        end
      end
  end
  return true
end

#TODO: rename get_non_interacting_orbs
" Return a set of all orbitals in os not affected by any kink. "
function get_non_interacting_orbs_of_set(Configuration::Configuration, os::Set{<:Orbital})
  non_int_orbs = Set{basis(Configuration)}()# use explicit constructor
  for orb in os
    if is_non_interacting(Configuration, orb)
      push!(non_int_orbs, orb)
    end
  end
  return non_int_orbs
end

# TODO: rename get_non_interacting_orbs_in_interval
" Return a set of all orbitals in os not affected by any kink in the open interval (τ_first,τ_last). "
function get_non_interacting_orbs_of_set_in_interval(Configuration::Configuration, os::Set{<:Orbital}, τ_first::ImgTime, τ_last::ImgTime )# :: Set{<:Orbital}
  non_int_orbs = Set{basis(Configuration)}()
  for orb in os
    if is_non_interacting_in_interval(Configuration, orb, τ_first, τ_last)
      push!(non_int_orbs, orb)
    end
  end
  return non_int_orbs
end


## Method definitions for function drop.
# This function is mostly used for calculating the changes proposed by an update
# in order to calculate proposal and/or acceptance probabilities.

" Return the configuration c without dropping anything. "
drop(c::Configuration, n::Nothing) = c

" Return a set with the orbital o dropped from oc. "
drop(oc::Set{T}, o::T) where {T <: Orbital} = setdiff(oc, Set([o]))
" Return a configuration with the orbital o dropped from c.occupations. "
drop(c::Configuration{T}, o::T) where {T <: Orbital} = Configuration(drop(c.occupations, Set([o])),c.kinks)

" Return a set with the orbitals in oc2 dropped from oc1. "
drop(oc1::Set{T}, oc2::Set{T}) where {T <: Orbital} = setdiff(oc1, oc2)
" Return a configuration with the orbitals in oc dropped from c.occupations. "
drop(c::Configuration{T}, oc::Set{T}) where {T <: Orbital} = Configuration(drop(c.occupations, oc), c.kinks)

" Return a SortedDict{ImgTime,<:Kink{<:Orbital}} with the pairs in ck2 dropped from ck1. "
drop(ck1::SortedDict{ImgTime,<:Kink{T}}, ck2::SortedDict{ImgTime,<:Kink{T}}) where {T <: Orbital} = setdiff(ck1, ck2)
" Return a configuration with the pairs in ck dropped from c.kinks. "
drop(c::Configuration{T}, ck::SortedDict{ImgTime,<:Kink{T}}) where {T <: Orbital} = Configuration(c.occupations, drop(c.kinks, ck))

" Return a SortedDict{ImgTime,<:Kink{<:Orbital}} with the pairs ps... dropped from ck. "
drop(ck::SortedDict{ImgTime,<:Kink{T}}, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital} = setdiff(ck, SortedDict(ps...))
" Return a configuration with the pairs ps... dropped from c.kinks. "
drop(c::Configuration{T}, ps::Pair{ImgTime,<:Kink{T}}...) where {T <: Orbital} = Configuration(c.occupations, drop(c.kinks, ps...))

" Return a configuration with occupations and kinks in c2 dropped from c1. "
drop(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital} = Configuration(drop(c1.occupations,c2.occupations), drop(c1.kinks,c2.kinks))


## Method definitions for function drop!.
# This function is mostly for applying the changes determined in an MC Step Δ
# to the current configuration c in function promote!(c, Δ). cf. CPIMC.jl
# These methods change the first argument in-place.

" Drop nothing. "
function drop!(c::Configuration, n::Nothing)
  nothing
end

" Drop a single orbital. "
function drop!(c::Configuration{T}, o::T) where {T <: Orbital}
    delete!(c.occupations, o)
end

" Drop all orbitals given in a set oc from c.occupations. "
function drop!(c::Configuration{T}, oc::Set{T}) where {T <: Orbital}
  for o in oc
    drop!(c, o)
  end
end

" Drop a single kink at time τ from c.kinks. "
function drop!(c::Configuration, τ::ImgTime)
  delete!(c.kinks, τ)
end

" Drop all kinks given in a SortedDict from c.kinks. "
function drop!(c::Configuration{T}, ck::SortedDict{ImgTime,<:Kink{T}}) where {T <: Orbital}
  for (τ, k) in ck
    drop!(c, τ)
  end
end

" Drop all kinks given by Varargs{Pair{ImgTime,<:Kink{T}}} from c.kinks. "
function drop!(c::Configuration{T}, pk::Pair{ImgTime,<:Kink{T}}...) where {T <: Orbital}
  for (τ, k) in pk
    drop!(c, τ)
  end
end

" Drop occupations and kinks given by second argument c2 from first argument c1. "
function drop!(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital}
  drop!(c1, c2.occupations)
  drop!(c1, c2.kinks)
end


## Method definitions for function add.
# This function is mostly used for calculating the changes proposed by an update
# in order to calculate proposal and/or acceptance probabilities.

" Return the configuration c without adding anything. "
add(c::Configuration, n::Nothing) = c

" Return a set with the orbital o added to oc. "
add(oc::Set{T}, o::T) where {T <: Orbital} = union(oc, Set([o]))
" Return a configuration with the orbital o added to c.occupations. "
add(c::Configuration{T}, o::T) where {T <: Orbital} = Configuration(add(c.occupations, o), c.kinks)

" Return a set with the orbitals in oc2 added to oc1. "
add(oc1::Set{T}, oc2::Set{T}) where {T <: Orbital} = union(oc1, oc2)
" Return a configuration with the orbitals in oc added to c.occupations. "
add(c::Configuration{T}, oc::Set{T}) where {T <: Orbital} = Configuration(union(c.occupations, oc), c.kinks)

" Return a SortedDict{ImgTime,<:Kink{<:Orbital}} with the pairs ps... added to ck. "
add(ck::SortedDict{ImgTime,<:Kink{T}}, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital} = merge(ck, SortedDict(ps...))
" Return a configuration with the pairs ps... added to c.kinks. "
add(c::Configuration{T}, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital} = Configuration(c.occupations, add(c.kinks, ps...))
" Return a configuration with the pairs in ck added to c.kinks. "
add(c::Configuration{T}, ck::SortedDict{ImgTime,<:Kink{T}}) where {T <: Orbital} = add(c, ck...)

" Return a configuration with occupations and kinks in c2 dropped from c1. "
add(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital} = Configuration(add(c1.occupations, c2.occupations), merge(c1.kinks, c2.kinks))


## Method definitions for function add!.
# This function is mostly for applying the changes determined in an MC Step Δ
# to the current configuration c in function promote!(c, Δ). cf. CPIMC.jl
# These methods change the first argument in-place.

" Add nothing. "
function add!(c::Configuration, n::Nothing)
  nothing
end

" Add a single orbital. "
function add!(c::Configuration{T}, o::T) where {T <: Orbital}
  push!(c.occupations, o)
end

" Add all orbitals given in a set oc to c.occupations. "
function add!(c::Configuration{T}, oc::Set{T}) where {T <: Orbital}
  union!(c.occupations, oc)
end

" Add a single kink at time τ to c.kinks. "
function add!(c::Configuration{T}, τ::ImgTime, k::Kink{T}) where {T <: Orbital}
  insert!(c.kinks, τ, k)
end

" Add all kinks given by the pairs ps... to ck::SortedDict{ImgTime,<:Kink{<:Orbital}}. "
add!(ck::SortedDict{ImgTime,<:Kink{T}}, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital} = merge!(ck, SortedDict(ps...))
" Add all kinks given by the pairs ps... to c.kinks. "
add!(c::Configuration{T}, ps::Pair{ImgTime,<:Kink{T}}...) where {T <: Orbital} = merge!(c.kinks, SortedDict(ps))

" Add all kinks given in a SortedDict to c.kinks. "
function add!(c::Configuration{T}, ck::SortedDict{ImgTime,<:Kink{T}}) where {T <: Orbital}
  merge!(c.kinks, ck)
end

" Add occupations and kinks given by second argument c2 from first argument c1. "
function add!(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital}
  add!(c1, c2.occupations)
  add!(c1, c2.kinks)
end







"""Returns True if left_kink and right_kink are entangled in a Type-B way.
This does not check wether the two kinks are neighbouring"""
function is_type_B(left_kink::T4, right_kink::T4)
  if (Set([left_kink.i, left_kink.j]) == Set([right_kink.k, right_kink.l])) &
        (Set([left_kink.k, left_kink.l]) == Set([right_kink.i, right_kink.j]))
    return(true)
  else
    return(false)
  end
end


"""Return a Tuple of 2 imaginaty times of "neighbouring" Kinks that are Type-B-Entangeld.
This function has no use in the current update set removability will therefore not be considered here.
"neighbouring" refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first. (for type-B-entanglement that always imples the
vice versa case)
The Set consists of the pairs where the Type-B-entanglement is oriented
to the left of the first τ."""
function get_left_type_B_pairs(c::Configuration)
  pairs_left = Set{Tuple{Fixed{Int64,60},Fixed{Int64,60}}}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = get_τ_borders(c, kink_orb_set ,τ)
    if is_type_B(c.kinks[τ_left], kink)
      push!(pairs_left, (τ, τ_left))
    end
  end
  return pairs_left
end


"""Return a Tuple of 2 imaginaty times of "neighbouring" Kinks that are Type-B-Entangeld. Removablility will
no be looked at. "neighbouring" refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first. (for type-B-entanglement that always imples the
vice versa case)
The Set consists of the pairs where the Type-B-entanglement is oriented
#to the right of the first τ."""
function get_right_type_B_pairs(c::Configuration)
  pairs_right = Set{Tuple{Fixed{Int64,60},Fixed{Int64,60}}}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = get_τ_borders(c, kink_orb_set ,τ)
    if is_type_B(kink, c.kinks[τ_right])
        push!(pairs_right, (τ, τ_right))
    end
  end
  return pairs_right
end


"""Return a Tuple of 2 imaginaty times of "neighbouring" Kinks that are Type-B-Entangeld AND removable.
"neighbouring" refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first. (for type-B-entanglement that always imples the
vice versa case)
The Set consists of the pairs where the Type-B-entanglement is oriented
#to the right of the first τ."""
function get_right_type_B_removable_pairs(c::Configuration)
  pairs_right = Set{Tuple{Fixed{Int64,60},Fixed{Int64,60}}}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = get_τ_borders(c, kink_orb_set ,τ)
    if is_type_B(kink, c.kinks[τ_right])
      if dot(kink.i.vec-kink.k.vec, kink.i.vec-kink.k.vec) <= ex_radius^2
        if kink.i.spin == kink.k.spin
          push!(pairs_right, (τ, τ_right))
        end
      end
    end
  end
  return pairs_right
end


"""Returns True if left_kink and right_kink are entangled in a Type-C way.
This does not check wether the two kinks are neighbouring"""
function is_type_C(left_kink::T4, right_kink::T4)
  if (Set([left_kink.i, left_kink.j]) == Set([right_kink.k, right_kink.l])) &
        !(Set([left_kink.k, left_kink.l]) == Set([right_kink.i, right_kink.j]))
    return(true)
  else
    return(false)
  end
end

"""Return a Tuple of 2 imaginaty times of "neighbouring" Kinks that are Type-C-Entangeld AND removable.
"neighbouring" refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-C-entanglement is oriented
#to the left of the first τ."""
function get_left_type_C_removable_pairs(c::Configuration)
  pairs_left = Set{Tuple{Fixed{Int64,60},Fixed{Int64,60}}}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = get_τ_borders(c, kink_orb_set ,τ)
    if is_type_C(c.kinks[τ_left], kink)
      if dot(kink.i.vec-kink.k.vec,kink.i.vec-kink.k.vec) <= (ex_radius^2)
        if kink.i.spin == kink.k.spin
          push!(pairs_left, (τ, τ_left))
        end
      end
    end
  end
  return pairs_left
end


"""Return a Tuple of 2 imaginaty times of "neighbouring" Kinks that are Type-C-Entangeld AND removable.
"neighbouring" refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-C-entanglement is oriented
#to the right of the first τ."""
function get_right_type_C_removable_pairs(c::Configuration)
  pairs_right = Set{Tuple{Fixed{Int64,60},Fixed{Int64,60}}}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = get_τ_borders(c, kink_orb_set ,τ)
    if is_type_C(kink, c.kinks[τ_right])
      if dot(kink.i.vec-kink.k.vec,kink.i.vec-kink.k.vec) <= (ex_radius^2)
        if kink.i.spin == kink.k.spin
          push!(pairs_right, (τ, τ_right))
        end
      end
    end
  end
  return pairs_right
end


"""Returns True if left_kink and right_kink are entangled in a Type-D way
This does not check wether the two kinks are neighbouring"""
function is_type_D(left_kink::T4, right_kink::T4)
  if !(Set([left_kink.i, left_kink.j]) == Set([right_kink.k, right_kink.l])) &
        (Set([left_kink.k, left_kink.l]) == Set([right_kink.i, right_kink.j]))
    return(true)
  else
    return(false)
  end
end

"""Return a Tuple of 2 imaginaty times of "neighbouring" Kinks that are Type-D-Entangeld AND removable.
"neighbouring" refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-D-entanglement is oriented
#to the left of the first τ."""
function get_left_type_D_removable_pairs(c::Configuration)
  pairs_left = Set{Tuple{Fixed{Int64,60},Fixed{Int64,60}}}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = get_τ_borders(c, kink_orb_set ,τ)
    if is_type_D(c.kinks[τ_left], kink)
      if dot(kink.i.vec-kink.k.vec,kink.i.vec-kink.k.vec) <= (ex_radius^2)
        if kink.i.spin == kink.k.spin
          push!(pairs_left, (τ, τ_left))
        end
      end
    end
  end
  return pairs_left
end

"""Return a Tuple of 2 imaginaty times of "neighbouring" Kinks that are Type-D-Entangeld AND removable.
"neighbouring" refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-D-entanglement is oriented
#to the right of the first τ."""
function get_right_type_D_removable_pairs(c::Configuration)
  pairs_right = Set{Tuple{Fixed{Int64,60},Fixed{Int64,60}}}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = get_τ_borders(c, kink_orb_set ,τ)
    if is_type_D(kink, c.kinks[τ_right])
      if dot(kink.i.vec-kink.k.vec,kink.i.vec-kink.k.vec) <= (ex_radius^2)
        if kink.i.spin == kink.k.spin
          push!(pairs_right, (τ, τ_right))
        end
      end
    end
  end
  return pairs_right
end


"""If left kink and right_kin are type-E-entangled it returns a tuple of the two kinks whos orbs
are sorted in a way that k an j of both kinks are the common orbitals.
Otherwise it returns false.
This does not check wether the two kinks are neighbouring"""
function is_type_E(left_kink::T4, right_kink::T4)
  c_orb1 = intersect!(Set([left_kink.i, left_kink.j]), Set([right_kink.k,right_kink.l]))
  c_orb2 = intersect!(Set([left_kink.k, left_kink.l]), Set([right_kink.i,right_kink.j]))
  if ((length(c_orb1) == 1) & (length(c_orb2) == 1))
    c_orb1 = first(c_orb1)
    c_orb2 = first(c_orb2)
    noncommon_orb_leftk_left = first(setdiff!(Set([left_kink.k, left_kink.l]), Set([right_kink.i,right_kink.j])))
    noncommon_orb_rightk_right = first(setdiff!(Set([right_kink.i,right_kink.j]), Set([left_kink.k, left_kink.l])))
    noncommon_orb_leftk_right = first(setdiff!(Set([left_kink.i, left_kink.j]), Set([right_kink.k,right_kink.l])))
    noncommon_orb_rightk_left = first(setdiff!(Set([right_kink.k,right_kink.l]), Set([left_kink.i, left_kink.j])))
    @assert(length(Set([noncommon_orb_leftk_left, noncommon_orb_rightk_right, noncommon_orb_leftk_right, noncommon_orb_rightk_left])) == 4)
    @assert(noncommon_orb_leftk_right.vec - noncommon_orb_leftk_left.vec == noncommon_orb_rightk_left.vec - noncommon_orb_rightk_right.vec)
    return((T4(noncommon_orb_leftk_right, c_orb1, c_orb2, noncommon_orb_leftk_left),
                T4(noncommon_orb_rightk_right, c_orb2, c_orb1, noncommon_orb_rightk_left)))
  else
    return(false)
  end
end

"""Return a Tuple of 2 imaginaty times of "neighbouring" Kinks that are Type-E-Entangeld AND removable.
"neighbouring" refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-E-entanglement is oriented
#to the left of the first τ."""
function get_left_type_E_removable_pairs(c::Configuration)
  pairs_left = Set{}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = get_τ_borders(c, kink_orb_set ,τ)
    sorted_kinks = is_type_E(c.kinks[τ_left], kink)
    if sorted_kinks == false
      continue
    end
    @assert(τ_left != τ_right)
    @assert(τ_right != ImgTime(1))
    if dot(last(sorted_kinks).i.vec - last(sorted_kinks).k.vec,
                last(sorted_kinks).i.vec - last(sorted_kinks).k.vec) <= (ex_radius^2)
      push!(pairs_left, (τ => last(sorted_kinks),τ_left => first(sorted_kinks)))
    end
  end
  return pairs_left
end


"""Return a Tuple of 2 imaginaty times of "neighbouring" Kinks that are Type-E-Entangeld AND removable.
"neighbouring" refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-E-entanglement is oriented
#to the right of the first τ."""
function get_right_type_E_removable_pairs(c::Configuration)
  pairs_right = Set{}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = get_τ_borders(c, kink_orb_set ,τ)

    sorted_kinks = is_type_E(kink, c.kinks[τ_right])
    if sorted_kinks == false
      continue
    end
    @assert(τ_left != τ_right)
    @assert(τ_right != ImgTime(1))
    if dot(first(sorted_kinks).i.vec - first(sorted_kinks).k.vec,
                first(sorted_kinks).i.vec - first(sorted_kinks).k.vec) <= (ex_radius^2)
        push!(pairs_right, (τ => first(sorted_kinks),τ_right => last(sorted_kinks)))
    end


  end
  return pairs_right
end


"""Returns an Array of the img-time ordered ladder operators used in the current configuration,
where an operator corresponds to a tuple of 1 (creator) or -1(annihilator) and the corresponding obrital"""
function get_timeordered_ladder_operators(c)
  operators = Array{Pair{Int,basis(c)},1}()
  for (_,kink) in c.kinks
    for orb in [kink.i,kink.j]
      push!(operators, 1 => orb)
    end
    for orb in [kink.k,kink.l]
      push!(operators, -1 => orb)
    end
  end
  return operators
end
"""Returns 1 or -1 depending on the order of all ladderoperators.
Used in the sign estimator."""
function get_ladder_operator_order_factor(c)
  phase_factor = 1
  operators = get_timeordered_ladder_operators(c)
  while !isempty(operators)
    if first(operators[1]) == 1
        op_type = 1
        index = 1
        op_orb = last(operators[index])

    else
        op_type = -1
        operators[1:2] = [operators[2], operators[1]]
        phase_factor *= -1
        index = 2
        op_orb = last(operators[index])
    end
    while operators[index + op_type] != (op_type*(-1) => op_orb)
      operators[index:index+1] = [operators[index+1], operators[index]]
      phase_factor *= -1
      index += 1
      @assert(index <= length(operators))
    end
    deleteat!(operators, sort([index + op_type,index]))
  end
  return (phase_factor)
end
