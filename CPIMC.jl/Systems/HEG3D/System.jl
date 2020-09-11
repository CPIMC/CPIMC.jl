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

"return quantum numbers of all states with energy ek"
function get_qstates(ek)
    # TODO: square-root of ek +1 ????
    d = Array{Array{Int,1},1}()
    kmax = floor(sqrt(ek))
    for x = -kmax:kmax
        for y = -kmax:kmax
            for z = -kmax:kmax
                if ek == x * x + y * y + z * z
                    push!(d, [x, y, z])
                end
            end
        end
    end
    d
end


function get_orblist_UEG(Emax)

    # estimate upper bound for maximum k-Component
    # kc = Int(ceil(sqrt(cutoff^(1.748626417827964) + 13*cutoff)))# cutoff radius (TODO: better estimate)

    # calculate energies of all possible k-vectors
    # where the absolute of each component is smaller than the cutoff value kc
    # the corresponding eigenenergies (without factor 1/2) are stored in container e
    K = DataFrame(x=Int64[], y=Int64[], z=Int64[], ek=Int64[], k=Array{Int64,1}[])
    K = DataFrame(k=vcat( collect( get_qstates(kk) for kk in 0:Emax) ) )
    filter!(row -> !isempty(row.k), K )# filter nonempty cols
    K.ek = dot(K.k,K.k)
    K.i = 1:size(K,1)
    return K
    # TODO: test sign implementation:
    # K = DataFrame(x=Int64[], y=Int64[], z=Int64[], ek=Int64[], k=Array{Int64,1}[], s=Sign[])
    ## assemble basis
    # for x in -kc:1:kc
    #     for y in -kc:1:kc
    #         for z in -kc:1:kc
    #             push!(K,[x,y,z,x*x+y*y+z*z, [x,y,z]])
    #             # TODO: test sign implementation:
    #             # push!(K,[x,y,z,x*x+y*y+z*z, [x,y,z]], Sign(false))
    #             # push!(K,[x,y,z,x*x+y*y+z*z, [x,y,z]], Sign(true))
    #         end
    #     end
    # end
    #
    # ## sort basis: lowest energies first
    # sort!(K,[:ek,:x,:y,:z])
    # # TODO: test sign implementation:
    # # sort!(K,[:ek,:x,:y,:z,:s])
    # ## add index column
    # K.i = 1:size(K,1)
    #
    # # return basis
    # K[K.i .<= cutoff,:]
end

function get_index(q::Array{Int64,1}, orblist::DataFrame)
    # orblist.i[(orblist.x .== q[1]) .& (orblist.y .== q[2]) .& (orblist.z .== q[3])][1]
    orblist.i[orblist.k .== [q]][1]
end

# function get_index(q::Array{Int64,1}, s::Sign, orblist::DataFrame)
    # orblist.i[(orblist.x .== q[1]) .& (orblist.y .== q[2]) .& (orblist.z .== q[3]) .& (orblist.s .== s)][1]
    # orblist.i[(orblist.k .== [q]) .& (orblist.s .== s)][1]
# end

function get_vector(index, orblist::DataFrame)
    orblist.k[orblist.i .== index][1]
end

function get_energy(index, orblist::DataFrame)
    orblist.ek[orblist.i .== index][1]
end

# function get_sign(index, orblist::DataFrame)
    # orblist.s[orblist.i .== index][1]
# end

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
