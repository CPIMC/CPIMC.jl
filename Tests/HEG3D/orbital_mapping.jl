
function get_index(q::Array{Int16,1}) :: Int64
    b = convert(Array{Int64,1}, q)
    b[1] << 32 | b[2] << 16 | b[3]
end

function get_vector(index) :: Array{Int16,1}
    convert(Array{Int16,1}, [ (index >>> 32) & 0xFFFF, (index >>> 16) & 0xFFFF, index & 0xFFFF ])
end

using DataFrames
# empty orblist
orblist = DataFrame(idx=Int64[], k=Array{Int16,1}[])

kmax = 20# maximum k component

for x in Int16.(-kmax:kmax)
    for y in Int16.(-kmax:kmax)
        for z in Int16.(-kmax:kmax)
            ik = Array{Int16,1}([x,y,z])
            ix = get_index(ik)
            push!(orblist, [ix, ik])#, get_vector(ix)])
        end
    end
end

println("Number of occurences of each index: ")
by(orblist, [:idx], count = :idx => length)

println("Number of occurences of each vector: ")
by(orblist, [:k], count = :k => length)

if any(by(orblist, [:idx], count = :idx => length).count .> 1)
    println("There are double indices!")
end
