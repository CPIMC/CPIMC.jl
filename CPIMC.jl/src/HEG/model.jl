using StaticArrays

import LinearAlgebra: dot

struct Ensemble
  "Brueckner parameter"
  rs :: Float64

  "reduced temperature"
  beta :: Float64

  "particle number"
  N :: Int
end

"representation of single particle state"
struct Orbital{D}
    "D-dimensional excitation vector"
    vec :: StaticVector{D,Int}
    "spin"
    spin :: Int
end

Orbital(v::Tuple,s=1) = Orbital(SVector(v),s)


function get_beta_internal(theta, N)
  return ((2*pi)^2)/(((6*(pi^2)*N)^(2/3))*theta)
end

function lambda(N::Int, rs::Float64)
    return (4/((2*pi)^3)) * (((4*pi)/3)^(1/3)) * rs * N^(1/3)
end

" single particle energy for momentum vector k "
function get_energy(o::Orbital)
    dot(o.vec,o.vec)
end



function get_abs_offdiagonal_element(e::Ensemble,c::Configuration,Kink::T4{Orbital{3}})
    @assert ((Kink.i.spin == Kink.k.spin) & (Kink.j.spin == Kink.l.spin))
    wijkl =  1/dot((Kink.i.vec-Kink.k.vec), (Kink.i.vec-Kink.k.vec))
    if Kink.i.spin == Kink.j.spin
        wijkl -= 1/dot((Kink.i.vec-Kink.l.vec), (Kink.i.vec-Kink.l.vec))
    end
    wijkl *= lambda(e.N,e.rs) * 1/2
    return abs(wijkl)
end

function get_orbshell(o::Orbital{1};dw::Int=2)
    eq = get_energy(o)
    qmax = Int(floor(eq))

    os = Set{Orbital{D}}()

    for x in -qmax:qmax
        if abs(x*x - eq) <= dw
            push!(os, Orbital((x,1)))
            push!(os, Orbital((x,-1)))
        end
    end

    os
end

function get_orbshell(o::Orbital{1};dw::Int=2)
    eq = get_energy(o)
    qmax = Int(floor(eq))

    os = Set{Orbital{2}}()

    for x in -qmax:qmax
        for y in -qmax:qmax
            if abs(x*x + y*y - eq) <= dw
                push!(os, Orbital((x,y),1))
                push!(os, Orbital((x,y),-1))
            end
        end
    end

    os
end

function get_orbshell(o::Orbital{3};dw::Int=2)
    eq = get_energy(o)
    qmax = Int(floor(eq))

    os = Set{Orbital{3}}()

    for x in -qmax:qmax
        for y in -qmax:qmax
            for z in -qmax:qmax
                if abs( x*x + y*y + z*z - eq ) <= dw
                    push!(os, Orbital((x,y,z),1))
                    push!(os, Orbital((x,y,z),-1))
                end
            end
        end
    end

    os
end

function get_sphere(o::Orbital{1}; dk::Int=2)
    os = Set{Orbital{1}}()

    for x in -dk:dk
        if x*x <= dk*dk
            push!(os, Orbital(o.vec+SVector(x),1))
            push!(os, Orbital(o.vec+SVector(x),-1))
        end
    end
    os
end

function get_sphere(o::Orbital{2}; dk::Int=2)
    os = Set{Orbital{2}}()

    for x in -dk:dk
        for y in -dk:dk
            if x*x + y*y <= dk*dk
                push!(os, Orbital(o.vec+SVector(x,y),1))
                push!(os, Orbital(o.vec+SVector(x,y),-1))
            end
        end
    end
    os
end


function get_sphere(o::Orbital{3}; dk::Int=1)
    os = Set{Orbital{3}}()

    for x in -dk:dk
        for y in -dk:dk
            for z in -dk:dk
                if x*x + y*y + z*z <= dk*dk
                    push!(os, Orbital(o.vec+SVector(x,y,z),1))
                    push!(os, Orbital(o.vec+SVector(x,y,z),-1))
                end
            end
        end
    end
    os
end

function get_orbs_with_spin(orbitals::Set{Orbital{3}},spin::Int)
    orbs_s = Set{Orbital{3}}()
    for orb in orbitals
        if orb.spin==spin
            push!(orbs_s, orb)
        end
    end
    return orbs_s
end

#returns a list of all Kinks that affect the given orbital
function get_kinks_of_orb(c::Configuration, orbital::Orbital)
  kinks_of_orb = SortedDict{img_time, Kink{Orbital{3}}}()
  for (tau_kink,kink) in c.kinks
    if (kink.i == orbital) | (kink.j == orbital) | (kink.k == orbital) | (kink.l == orbital)
      kinks_of_orb[tau_kink] = kink
    end
  end
  return kinks_of_orb
end

#returns a tuple of the imiginary times between 0 and 1 of the nearest Kinks before Tau and the nearest Kink after Tau
#has to be multiplied with beta to get reeal imaginary times
function get_nearest_Tau_effecting_orb(Configuration::Configuration, orbital::Orbital,Tau::img_time)
  current_tau = 0
  Kinks_of_orb = get_kinks_of_orb(Configuration, orbital)
  if length(Kinks_of_orb) == 0
      return(0,1)
  end
  #TO DO: Binäre Suche benutzten
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

function get_Tau_boarders(Configuration::Configuration, orbitals::Set{Orbital{3}},Tau::img_time)
  if length(Configuration.kinks) == 0
      return (0,1)
  end
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

  for orb in orbitals
    Tupel = get_nearest_Tau_effecting_orb(Configuration, orb, Tau)

    #hier muss man immer prüfen ob das Intervall die Grenze Beta bzw Null überschreitet
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
  return (Tau_left,Tau_right)
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

#returns all orbs with no kinks
function get_non_interacting_orbs_of_set(Configuration::Configuration, os::Set{Orbital{3}})
  non_int_orbs = Set{Orbital{3}}()
  for orb in os
    if is_non_interacting(Configuration, orb)
      push!(non_int_orbs, orb)
    end
  end
  return(non_int_orbs)
end
