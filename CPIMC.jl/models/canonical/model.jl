
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

function get_beta_internal(theta, N)
  return ((2*pi)^2)/(((6*(pi^2)*N)^(2/3))*theta)
end

mutable struct Configuration
  "set of currently occupied orbitals"
  occupations :: Set{Int}
end

function emptyOrbs(e::Ensemble, c::Configuration)
  return filter(x -> !(x in c.occupations), 1:e.cutoff)
end
