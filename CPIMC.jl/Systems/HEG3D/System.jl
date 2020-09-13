## Calculate a mapping from orbital indices to orbital quantum numbers
## e.g. kvectors
## The units are so that the possible values of each k-vector component are integer numbers (i.e. divide by 2pi/boxarea compared to a.u.)

using DataFrames
import LinearAlgebra: dot

# wrapper for the sign
# struct Sign
#     sign :: Bool
# end
# TODO : Avoid fields with abstract type
# mutable struct Sign{T<:Bool}
    # sign :: T
# end

function get_index(q::Array{Int16,1}) :: Int64
    b = convert(Array{Int64,1}, q)
    b[1] << 32 | b[2] << 16 | b[3]
end

function get_vector(index) :: Array{Int16,1}
    convert(Array{Int16,1}, [ (index >>> 32) & 0xFFFF, (index >>> 16) & 0xFFFF, index & 0xFFFF ])
end

function get_energy(index)
    vec = get_vector(index)
    dot(vec,vec)
end

### Estimators
function Ekin(e::Ensemble, c::Configuration, orblist::DataFrame)
    sum(get_energy(n,orblist) for n in c.occupations)
end

function occVec(e::Ensemble, c::Configuration, orblist::DataFrame)
  return map(x -> Int(x in c.occupations), 1:e.cutoff)
end

function particleNumber(e::Ensemble, c::Configuration, orblist::DataFrame)
  return c.N
end


### Units
function get_beta_internal(theta, N)
  return ((2*pi)^2)/(((6*(pi^2)*N)^(2/3))*theta)
end
