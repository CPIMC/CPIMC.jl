using DataStructures


struct T2{T}
  " start orbital "
  i :: T

  " end orbital "
  j :: T
end

struct T4{T}
  " start orbitals "
  i :: T
  j :: T

  " end orbitals "
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


