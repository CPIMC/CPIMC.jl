## Calculate a mapping from orbital indices to orbital quantum numbers
## e.g. kvectors
## The units are so that the possible values of each k-vector component are integer numbers (i.e. divide by 2pi/boxarea compared to a.u.)

using DataFrames

# wrapper for the sign
# struct Sign
#     sign :: Bool
# end
# TODO : Avoid fields with abstract type
# mutable struct Sign{T<:Bool}
    # sign :: T
# end

function get_orblist_UEG(cutoff)

    # estimate upper bound for maximum k-Component
    kc = Int(ceil(sqrt(cutoff^(1.748626417827964) + 13*cutoff)))# cutoff radius (TODO: better estimate)

    # calculate energies of all possible k-vectors
    # where the absolute of each component is smaller than the cutoff value kc
    # the corresponding eigenenergies (without factor 1/2) are stored in container e
    K = DataFrame(x=Int64[], y=Int64[], z=Int64[], ek=Int64[], k=Array{Int64,1}[])
    # TODO: test sign implementation:
    # K = DataFrame(x=Int64[], y=Int64[], z=Int64[], ek=Int64[], k=Array{Int64,1}[], s=Sign[])
    ## assemble basis
    for x in -kc:1:kc
        for y in -kc:1:kc
            for z in -kc:1:kc
                push!(K,[x,y,z,x*x+y*y+z*z, [x,y,z]])
                # TODO: test sign implementation:
                # push!(K,[x,y,z,x*x+y*y+z*z, [x,y,z]], Sign(false))
                # push!(K,[x,y,z,x*x+y*y+z*z, [x,y,z]], Sign(true))
            end
        end
    end

    ## sort basis: lowest energies first
    sort!(K,[:ek,:x,:y,:z])
    # TODO: test sign implementation:
    # sort!(K,[:ek,:x,:y,:z,:s])
    ## add index column
    K.i = 1:size(K,1)

    # return basis
    K[K.i .<= cutoff,:]
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

# function get_sign(index, orblist::DataFrame)
    # orblist.s[orblist.i .== index][1]
# end
