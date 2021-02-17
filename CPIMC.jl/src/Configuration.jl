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
      return ("nix","nix")
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
  # Initially we set τ right and τl left to the nearest kinks left and right of τ.
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
    if tupel[1] == "nix"
        non_interacting_orb_counter += 1
    else
        # here we always have to check wether the given intervall extends over 1
        if τ_left < τ < tupel[1]
            # "nix"
        elseif tupel[1] < τ < τ_left
            τ_left = tupel[1]
        elseif τ_left < tupel[1]
          τ_left = tupel[1]
        end

        if tupel[2] < τ < τ_right
            # "nix"
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
            # "nix"
        elseif (kink.i == orbital) | (kink.j == orbital) | (kink.k == orbital) | (kink.l == orbital)
              return false
        end
      end
  else
      for (τ_kink,kink) in Configuration.kinks
        if ((τ_kink <= τ_first) & (τ_kink >= τ_last))
            # "nix"
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


## Method definitions for function add.
# This function is mostly used for calculating the changes proposed by an update
# in order to calculate proposal and/or acceptance probabilities.

" Return the configuration c without dropping anything. "
drop(c::Configuration, n::Nothing) = c

" Drop a single orbital. "
function drop!(c::Configuration{T}, o::T) where {T <: Orbital}
  delete!(c.occupations, o)
end

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
  # for k in ck
  #   add!(c, k...)
  # end
  merge!(c.kinks, ck)
end

" Add occupations and kinks given by second argument c2 from first argument c1. "
function add!(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital}
  add!(c1, c2.occupations)
  add!(c1, c2.kinks)
end
