
using StaticArrays

struct Ensemble
  "number of basis orbitals"
  cutoff :: Int

  "Brueckner parameter"
  rs :: Float64

  "reduced temperature"
  beta :: Float64

  "particle number"
  N :: Int
end

struct Orbital{T}
    "single particle quantum numbers of spinless ueg"
    qnums :: StaticVector{T,Int16}
end

# struct Basis
#   "set of all orbitals in the basis"
#   Orbitals :: Set{Orbital}
# end

mutable struct Configuration
  "set of currently occupied orbitals"
  occupations :: Set{Orbital}
end

function emptyOrbs(c::Configuration, o::Set{Orbital}) :: Set{Orbital}
    setdiff!(c.occupations, o)
end
