using DataStructures


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


" multi-particle trajectory using single particle states with type T "
mutable struct Configuration{T}
  " set of orbitals occupied at tau=0 "
  occupations :: Set{T}

  " excitations, using tau as an index "
  kinks :: SortedDict{Float64, Kink{T}}

  Configuration(s::Set{T}) where {T} = new{T}(s,SortedDict{Float64,Kink}(Base.Forward))
end

function change_occupations(occs::Set, K::T4)
  assert(in(K.k, occs) & in(K.l, occs))
  assert(!in(K.i, occs) & !in(K.j, occs))
  delete!(occs, k)
  delete!(occs, l)
  push!(occs, i)
  push!(occs, j)
end

function get_occupations_at(conf::Configuration, Tau::Float64)
  occupations = copy(conf.occupations)
  for (tau_kink,kink) in conf.kinks
    if tau_kink < Tau
      change_occupations(occupations, kink)
    else
      break
    end
  end
  return occupations
end

function get_kinks_of(Configuration::Configuration, orbital::Orbital)
  kinks = SortedDict{}
  for (tau_kink,kink) in Configuration.kinks
    if kink.i == orbital | kink.j == orbital | kink.k == orbital | kink.l == orbital
      kinks{Tau_Kink} = kink
    end
  end
  return kinks
end

function is_non_interacting(Configuration::Configuration, orbital::Orbital)
  for tau_kink,kink in Configuration.kinks
    if kink.i == orbital | kink.j == orbital | kink.k == orbital | kink.l == orbital
      return(false)
    end
  end
  return(true)
end


function get_non_interacting_orbs_of_set(Configuration, os::Set{Orbital})
  non_int_orbs = Set{Orbital}()
  for orb in os:
    if is_non_interacting(Configuration, orb)
      push!(non_int_orbs, orb)
    end
  end
  return(non_int_orbs)
end
