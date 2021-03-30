abstract type Ensemble end

struct CEnsemble <: Ensemble
  "coupling parameter"
  λ :: Float64
  "inverse temperature"
  β :: Float64
  "particle number"
  N :: Int
end


struct GCEnsemble <: Ensemble
  "coupling parameter"
  λ :: Float64
  "inverse temperature"
  β :: Float64
  "chemical potential"
  µ :: Float64
end


