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
  Configuration(s::Set{T}, k:: SortedDict{Float64, Kink{T}}) where {T} = new{T}(s,k)
end

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
