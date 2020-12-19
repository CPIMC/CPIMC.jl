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


const img_time = Fixed{Int64,60}

" multi-particle trajectory using single particle states with type T "
mutable struct Configuration{T}
  " set of orbitals occupied at tau=0 "
  occupations :: Set{T}

  " excitations, using tau as an index "
  kinks :: SortedDict{img_time, Kink{T}}

  #constructor für eine Ideale konfiguration
  Configuration(s::Set{T}) where {T} = new{T}(s,SortedDict{img_time,Kink}(Base.Forward))
  Configuration(s::Set{T}, k::SortedDict{img_time, Kink{T}}) where {T} = new{T}(s,k)
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

#Find occupation numbers at Tau if there is a Kink at Tau find occupations right from it
function get_occupations_at(conf::Configuration, Tau::img_time)
  occupations = copy(conf.occupations)
  for (tau_kink,kink) in conf.kinks
    if tau_kink <= Tau
      change_occupations(occupations, kink)
    else
      break
    end
  end
  return occupations
end


#returns a list of all Kinks that affect the given orbital
function get_kinks_of_orb(c::Configuration, orbital::Orbital)
  kinks_of_orb = SortedDict{img_time, Kink{<:Orbital}}()
  for (tau_kink,kink) in c.kinks
    if (kink.i == orbital) | (kink.j == orbital) | (kink.k == orbital) | (kink.l == orbital)
      kinks_of_orb[tau_kink] = kink
    end
  end
  return kinks_of_orb
end

#Returns a tuple of the imiginary times between 0 and 1 of the nearest Kinks before Tau and the nearest Kink after Tau,
# which effect the given orbital. Has to be multiplied with beta to get real imaginary times.
function get_nearest_Tau_effecting_orb(Configuration::Configuration, orbital::Orbital,Tau::img_time)
  current_tau = 0
  Kinks_of_orb = get_kinks_of_orb(Configuration, orbital)
  if length(Kinks_of_orb) == 0
      return("nix","nix")
  end
  #TO DO: Binäre Suche benutzten?
  for (tau_kink,kink) in Kinks_of_orb
      if tau_kink > Tau
        if current_tau == 0
          return (first(last(Kinks_of_orb)),tau_kink)
        else
          return (current_tau,tau_kink)
        end
      else
          #es soll nicht Tau als Grenze zurückgegeben werden
          if tau_kink != Tau
              current_tau = tau_kink
          end
      end
  end
  return (current_tau, first(first(Kinks_of_orb)))
end



function get_Tau_boarders(Configuration::Configuration, orbitals::Set{<:Orbital},Tau::img_time)
  if length(Configuration.kinks) == 0
      return(img_time(0),img_time(1))
  end
  #Initially we set Tau right and Taul left to the nearest Kinks left and right of Tau.
  Tau_left_semi_token  = searchsortedafter(Configuration.kinks, Tau)
  Tau_right_semi_token = searchsortedlast(Configuration.kinks, Tau)
  if Tau_left_semi_token == pastendsemitoken(Configuration.kinks)
      Tau_left = first(first(Configuration.kinks))
  else
      Tau_left = first(deref((Configuration.kinks, Tau_left_semi_token)))
  end
  if Tau_right_semi_token == beforestartsemitoken(Configuration.kinks)
      Tau_right = first(last(Configuration.kinks))
  else
      Tau_right = first(deref((Configuration.kinks, Tau_right_semi_token)))
      if Tau_right == Tau
          Tau_right_semi_token = regress((Configuration.kinks,Tau_right_semi_token))
          if Tau_right_semi_token == beforestartsemitoken(Configuration.kinks)
              Tau_right = first(last(Configuration.kinks))
          else
              Tau_right = first(deref((Configuration.kinks, Tau_right_semi_token)))
          end
      end
  end
  #Now search for the nearst Taus that do actually effect one of the Orbitals
  non_interacting_orb_counter = 0
  for orb in orbitals
    @assert(Tau_right != img_time(1))
    Tupel = get_nearest_Tau_effecting_orb(Configuration, orb, Tau)
    if Tupel[1] == "nix"
        non_interacting_orb_counter += 1
    else
        #here we always have to check wether the given intervall extends over 1
        if Tau_left < Tau < Tupel[1]
            "nix"
        elseif Tupel[1] < Tau < Tau_left
            Tau_left = Tupel[1]
        elseif Tau_left < Tupel[1]
          Tau_left = Tupel[1]
        end

        if Tupel[2] < Tau < Tau_right
            "nix"
        elseif Tau_right < Tau < Tupel[2]
          Tau_right = Tupel[2]
        elseif Tupel[2] < Tau_right
            Tau_right = Tupel[2]
        end
    end
  end
  if non_interacting_orb_counter == length(orbitals)
      return return(img_time(0),img_time(1))
  else
      return (Tau_left,Tau_right)
  end
end



#see if an orb has no Kinks
function is_non_interacting(Configuration::Configuration, orbital::Orbital)
  for (tau_kink,kink) in Configuration.kinks
    if (kink.i == orbital) | (kink.j == orbital) | (kink.k == orbital) | (kink.l == orbital)
      return(false)
    end
  end
  return(true)
end

#see if an orb has no Kinks between two Taus (ignoring Kinks at one of the Taus)
function is_non_interacting_in_interval(Configuration::Configuration, orbital::Orbital, Tau_first::img_time, Tau_last::img_time)
  @assert Tau_first != Tau_last
  if Tau_first < Tau_last
      for (tau_kink,kink) in Configuration.kinks
        if (tau_kink <= Tau_first) | (tau_kink >= Tau_last)
            "nix"
        elseif (kink.i == orbital) | (kink.j == orbital) | (kink.k == orbital) | (kink.l == orbital)
              return(false)
        end
      end
  else
      for (tau_kink,kink) in Configuration.kinks
        if ((tau_kink <= Tau_first) & (tau_kink >= Tau_last))
            "nix"
        elseif (kink.i == orbital) | (kink.j == orbital) | (kink.k == orbital) | (kink.l == orbital)
              return(false)
        end
      end
  end
  return(true)
end

#returns all orbs with no kinks
function get_non_interacting_orbs_of_set(Configuration::Configuration, os::Set{<:Orbital})
  non_int_orbs = Set() #Set{<:Orbital}() funktioniert nicht, anscheineinend lassen
  #sich Set-objekte nicht mit angabe eine abstarken types erstellen.
  for orb in os
    if is_non_interacting(Configuration, orb)
      push!(non_int_orbs, orb)
    end
  end
  return(non_int_orbs)
end

#returns all orbs with no kinks between 2 taus, ignoring Kinks at one of the Taus
function get_non_interacting_orbs_of_set_in_interval(Configuration::Configuration, os::Set{<:Orbital}, Tau_first::img_time, Tau_last::img_time )
  non_int_orbs = Set() #Set{<:Orbital}() funktioniert nicht, anscheineinend lassen
  #sich Set-objekte nicht mit angabe eine abstarken types erstellen.
  for orb in os
    if is_non_interacting_in_interval(Configuration, orb, Tau_first, Tau_last)
      push!(non_int_orbs, orb)
    end
  end
  return(non_int_orbs)
end
