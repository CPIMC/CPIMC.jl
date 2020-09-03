#Das noch großkanonisch
struct Ensemble
  "number of basis orbitals"
  cutoff :: Int

  "box length"
  rs :: Float64

  "redused temperature"
  theta :: Float64

  "particle number"
  N :: Int
end

mutable struct Configuration
  "set of currently occupied orbitals"
  occupations :: Set{Int}

  "number of particles"
  N :: Int
end

function K(n::Int, e::Ensemble, c::Configuration)
  return pi^2 / (2 * e.L^2) * n^2
end

function emptyOrbs(e::Ensemble, c::Configuration)
  return filter(x -> !(x in c.occupations), 1:e.cutoff)
end

function occVec(e::Ensemble, c::Configuration)
  return map(x -> Int(x in c.occupations), 1:e.cutoff)
end

function totalEnergy(e::Ensemble, c::Configuration)
  return sum(map(x -> K(x,e,c), collect(c.occupations)))
end

function particleNumber(e::Ensemble, c::Configuration)
  return c.N
end
