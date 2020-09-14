
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

mutable struct Configuration
  "set of currently occupied orbitals"
  occupations :: Set{Int}
end
