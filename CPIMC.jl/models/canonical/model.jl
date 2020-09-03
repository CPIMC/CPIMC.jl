
struct Ensemble
  "number of basis orbitals"
  cutoff :: Int

  "Brueckner parameter"
  rs :: Float64

  "reduced temperature"
  theta :: Float64

  "particle number"
  N :: Int
end

function get_beta_internal(theta, N)
  return ((2*pi)^2)/(((6*(pi^2)*N)^(2/3))*theta)
end

mutable struct Configuration
  "set of currently occupied orbitals"
  occupations :: Set{Int}
end
