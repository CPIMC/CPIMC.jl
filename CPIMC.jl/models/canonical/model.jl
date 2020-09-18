
using StaticArrays

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

Orbital(v::Tuple,s=0) = Orbital(SVector(v),s)

mutable struct Configuration
  "set of currently occupied orbitals"
  occupations :: Set{Orbital}
end
