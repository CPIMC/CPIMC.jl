
function get_index(q::Array{Int16,1}) :: Int64
    b = convert(Array{Int64,1}, q)
    b[1] << 32 | b[2] << 16 | b[3]
end

function get_vector(index) :: Array{Int16,1}
    convert(Array{Int16,1}, [ (index >>> 32) & 0xFFFF, (index >>> 16) & 0xFFFF, index & 0xFFFF ])
end
using DataFrames

orblist = DataFrame(idx=Int64[], k=Array{Int16,1}[])#, calc_k=Array{Int16,1}[])

kmax = 20

for x in Int16.(-kmax:kmax)
    for y in Int16.(-kmax:kmax)
        for z in Int16.(-kmax:kmax)
            ik = Array{Int16,1}([x,y,z])
            ix = get_index(ik)
            push!(orblist, [ix, ik])#, get_vector(ix)])
        end
    end
end


orblist
