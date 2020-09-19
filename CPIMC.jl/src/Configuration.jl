" multi-particle trajectory using single particle states with type T "
mutable struct Configuration{T}
  " set of orbitals occupied at tau=0 "
  occupations :: Set{T}
end
