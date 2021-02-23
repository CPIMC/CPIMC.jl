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


const ImgTime = Fixed{Int64,60}

" multi-particle trajectory using single particle states with type T "
mutable struct Configuration{T}
  " set of orbitals occupied at τ=0 "
  occupations :: Set{T}
  " excitations, using τ as an index "
  kinks :: SortedDict{ImgTime, Kink{T}}

  sign :: Int8

  #constructor für eine Ideale konfiguration
  Configuration(o::Set{T}) where {T} = new{T}(o,SortedDict{ImgTime,Kink}(Base.Forward),1)
  Configuration(o::Set{T}, k::SortedDict{ImgTime, Kink{T}}) where {T} = new{T}(o,k,1)
  Configuration(o::Set{T}, k::SortedDict{ImgTime, Kink{T}}, s::Int8) where {T} = new{T}(o,k,s)
end

abstract type Orbital end


#Execute a type-4 Kink on a set of occupationnumbers
function change_occupations(occs::Set, K::T4)
  #try
  @assert (in(K.k, occs) & in(K.l, occs))
  @assert (!in(K.i, occs) & !in(K.j, occs))
  #catch ex
  #    print("BREAK")
  #end
  delete!(occs, K.k)
  delete!(occs, K.l)
  push!(occs, K.i)
  push!(occs, K.j)
end

#Find occupation numbers at τ if there is a Kink at τ find occupations right from it
function get_occupations_at(conf::Configuration, τ::ImgTime)
  occupations = copy(conf.occupations)
  for (τ_kink,kink) in conf.kinks
    if τ_kink <= τ
      change_occupations(occupations, kink)
    else
      break
    end
  end
  return occupations
end


#returns a list of all Kinks that affect the given orbital
function get_kinks_of_orb(c::Configuration, orbital::Orbital)
  kinks_of_orb = SortedDict{ImgTime, Kink{<:Orbital}}()
  for (τ_kink,kink) in c.kinks
    if (kink.i == orbital) | (kink.j == orbital) | (kink.k == orbital) | (kink.l == orbital)
      kinks_of_orb[τ_kink] = kink
    end
  end
  return kinks_of_orb
end

#Returns a tuple of the imaginary times between 0 and 1 of the nearest Kinks before τ and the nearest Kink after τ,
# which effect the given orbital. Has to be multiplied with β to get real imaginary times.
function get_nearest_τ_affecting_orb(Configuration::Configuration, orbital::Orbital,τ::ImgTime)
  current_τ = 0
  kinks_of_orb = get_kinks_of_orb(Configuration, orbital)
  if length(kinks_of_orb) == 0
      return("nix","nix")
  end
  #TO DO: Binäre Suche benutzten?
  for (τ_kink,kink) in kinks_of_orb
      if τ_kink > τ
        if current_τ == 0
          return (first(last(kinks_of_orb)),τ_kink)
        else
          return (current_τ,τ_kink)
        end
      else
          #es soll nicht τ als Grenze zurückgegeben werden
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
  #Now search for the nearst τs that do actually effect one of the Orbitals
  non_interacting_orb_counter = 0
  for orb in orbitals
    @assert(τ_right != ImgTime(1))
    tupel = get_nearest_τ_affecting_orb(Configuration, orb, τ)
    if tupel[1] == "nix"
        non_interacting_orb_counter += 1
    else
        #here we always have to check wether the given intervall extends over 1
        if τ_left < τ < tupel[1]
            "nix"
        elseif tupel[1] < τ < τ_left
            τ_left = tupel[1]
        elseif τ_left < tupel[1]
          τ_left = tupel[1]
        end

        if tupel[2] < τ < τ_right
            "nix"
        elseif τ_right < τ < tupel[2]
          τ_right = tupel[2]
        elseif tupel[2] < τ_right
            τ_right = tupel[2]
        end
    end
  end
  if non_interacting_orb_counter == length(orbitals)
      return return(ImgTime(0),ImgTime(1))
  else
      return (τ_left,τ_right)
  end
end



#see if an orb has no Kinks
function is_non_interacting(Configuration::Configuration, orbital::Orbital)
  for (τ_kink,kink) in Configuration.kinks
    if (kink.i == orbital) | (kink.j == orbital) | (kink.k == orbital) | (kink.l == orbital)
      return(false)
    end
  end
  return(true)
end

#see if an orb has no Kinks between two τs (ignoring Kinks at one of the τs)
function is_non_interacting_in_interval(Configuration::Configuration, orbital::Orbital, τ_first::ImgTime, τ_last::ImgTime)
  @assert τ_first != τ_last
  if τ_first < τ_last
      for (τ_kink,kink) in Configuration.kinks
        if (τ_kink <= τ_first) | (τ_kink >= τ_last)
            "nix"
        elseif (kink.i == orbital) | (kink.j == orbital) | (kink.k == orbital) | (kink.l == orbital)
              return(false)
        end
      end
  else
      for (τ_kink,kink) in Configuration.kinks
        if ((τ_kink <= τ_first) & (τ_kink >= τ_last))
            "nix"
        elseif (kink.i == orbital) | (kink.j == orbital) | (kink.k == orbital) | (kink.l == orbital)
              return(false)
        end
      end
  end
  return(true)
end

#returns all orbs with no kinks
function get_non_interacting_orbs_of_set(Configuration::Configuration, os::Set{<:Orbital})# ":: Set{<:Orbital}" funktioniert nicht
  non_int_orbs = Set() #Set{<:Orbital}() funktioniert nicht, anscheineinend lassen
  #sich Set-objekte nicht mit angabe eine abstarken types erstellen.
  for orb in os
    if is_non_interacting(Configuration, orb)
      push!(non_int_orbs, orb)
    end
  end
  return(non_int_orbs)
end

#returns all orbs with no kinks between 2 τs, ignoring Kinks at one of the τs
function get_non_interacting_orbs_of_set_in_interval(Configuration::Configuration, os::Set{<:Orbital}, τ_first::ImgTime, τ_last::ImgTime )# :: Set{<:Orbital}
  non_int_orbs = Set() #Set{<:Orbital}() funktioniert nicht, anscheineinend lassen
  #sich Set-objekte nicht mit angabe eine abstarken types erstellen.
  for orb in os
    if is_non_interacting_in_interval(Configuration, orb, τ_first, τ_last)
      push!(non_int_orbs, orb)
    end
  end
  return(non_int_orbs)
end

#Returns True if left_kink and right_kink are entangled in a Type-B way
function is_type_B(left_kink::T4, right_kink::T4)
  if (Set([left_kink.i, left_kink.j]) == Set([right_kink.k, right_kink.l])) &
        (Set([left_kink.k, left_kink.l]) == Set([right_kink.i, right_kink.j]))
    return(true)
  else
    return(false)
  end
end

#Returns True if left_kink and right_kink are entangled in a Type-C way
function is_type_C(left_kink::T4, right_kink::T4)
  if (Set([left_kink.i, left_kink.j]) == Set([right_kink.k, right_kink.l])) &
        !(Set([left_kink.k, left_kink.l]) == Set([right_kink.i, right_kink.j]))
    return(true)
  else
    return(false)
  end
end

#Returns True if left_kink and right_kink are entangled in a Type-D way
function is_type_D(left_kink::T4, right_kink::T4)
  if !(Set([left_kink.i, left_kink.j]) == Set([right_kink.k, right_kink.l])) &
        (Set([left_kink.k, left_kink.l]) == Set([right_kink.i, right_kink.j]))
    return(true)
  else
    return(false)
  end
end


#If left kink and right_kin are type-e-entangled it returns a tuple of the two kinks whos orbs
#are sorted in a way that k an j of both kinks are common orbitals
#otherwise returns false
#This does only work if at least one of the kinks is a neightboir of the other
function is_type_E(left_kink::T4, right_kink::T4)
  c_orb1 = intersect!(Set([left_kink.i, left_kink.j]), Set([right_kink.k,right_kink.l]))
  c_orb2 = intersect!(Set([left_kink.k, left_kink.l]), Set([right_kink.i,right_kink.j]))
  """if (length(intersect!(Set([left_kink.i, left_kink.j]), Set([right_kink.i,right_kink.j]))) != 0) |
        (length(intersect!(Set([left_kink.k, left_kink.l]), Set([right_kink.k,right_kink.l]))) != 0)
        return(false)
  end"""
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



#Return a two(a Tuple of) Sets of Tuples of "neighbouring" Kink that are Type-B-Entangeld.
#"neighbouring" refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
#other kink in the corresponding direktion are looked at.
#The Tuples are always arranged in a way that the Kink who gets neighboured by
#the opther stands first. (for type-B-entanglement that always imples the
#vice versa case)
#The Set consists of the pairs where the Type-B-entanglement is oriented
#to the left of the first τ.
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

#Return a two(a Tuple of) Sets of Tuples of "neighbouring" Kink that are Type-B-Entangeld.
#"neighbouring" refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
#other kink in the corresponding direktion are looked at.
#The Tuples are always arranged in a way that the Kink who gets neighboured by
#the opther stands first. (for type-B-entanglement that always imples the
#vice versa case)
#The Set consists of the pairs where the Type-B-entanglement is oriented
#to the right of the first τ.
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

#Return a two(a Tuple of) Sets of Tuples of "neighbouring" Kink that are Type-B-Entangeld.
#"neighbouring" refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
#other kink in the corresponding direktion are looked at.
#The Tuples are always arranged in a way that the Kink who gets neighboured by
#the opther stands first. (vice versa does not have to be the case)
#The Set consists of the pairs where the Type-C-entanglement is oriented
#to the left of the first τ.
function get_left_type_C_pairs(c::Configuration)
  pairs_left = Set{Tuple{Fixed{Int64,60},Fixed{Int64,60}}}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = get_τ_borders(c, kink_orb_set ,τ)
    if is_type_C(c.kinks[τ_left], kink)
      if dot(kink.i.vec-kink.k.vec,kink.i.vec-kink.k.vec) <= (ex_radius^2)
        push!(pairs_left, (τ, τ_left))
      end
    end
  end
  return pairs_left
end


#Return a two(a Tuple of) Sets of Tuples of "neighbouring" Kink that are Type-B-Entangeld.
#"neighbouring" refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
#other kink in the corresponding direktion are looked at.
#The Tuples are always arranged in a way that the Kink who gets neighboured by
#the opther stands first. (vice versa does not have to be the case)
#The Set consists of the pairs where the Type-C-entanglement is oriented
#to the right of the first τ.
function get_right_type_C_pairs(c::Configuration)
  pairs_right = Set{Tuple{Fixed{Int64,60},Fixed{Int64,60}}}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = get_τ_borders(c, kink_orb_set ,τ)
    if is_type_C(kink, c.kinks[τ_right])
      if dot(kink.i.vec-kink.k.vec,kink.i.vec-kink.k.vec) <= (ex_radius^2)
        push!(pairs_right, (τ, τ_right))
      end
    end
  end
  return pairs_right
end


function get_left_type_D_pairs(c::Configuration)
  pairs_left = Set{Tuple{Fixed{Int64,60},Fixed{Int64,60}}}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = get_τ_borders(c, kink_orb_set ,τ)
    if is_type_D(c.kinks[τ_left], kink)
      if dot(kink.i.vec-kink.k.vec,kink.i.vec-kink.k.vec) <= (ex_radius^2)
        push!(pairs_left, (τ, τ_left))
      end
    end
  end
  return pairs_left
end


function get_right_type_D_pairs(c::Configuration)
  pairs_right = Set{Tuple{Fixed{Int64,60},Fixed{Int64,60}}}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = get_τ_borders(c, kink_orb_set ,τ)
    if is_type_D(kink, c.kinks[τ_right])
      if dot(kink.i.vec-kink.k.vec,kink.i.vec-kink.k.vec) <= (ex_radius^2)
        push!(pairs_right, (τ, τ_right))
      end
    end
  end
  return pairs_right
end

#Get Kinks that are Type-E-removable and were the right kink is neighbouring the leftone
#the first partt of the returned tupel corresponds to the left kink
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

#Get Kinks that are Type-E-removable and were the left kink is neighbouring the right one
#the first partt of the returned tupel corresponds to the right kink
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
