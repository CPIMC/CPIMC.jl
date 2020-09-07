import LinearAlgebra: norm

function get_degen(ek, kmax)
    d = 0
    for x = -kmax:kmax
        for y = -kmax:kmax
            for z = -kmax:kmax
                if ek == x * x + y * y + z * z
                    d += 1
                end
            end
        end
    end
    d
end

function get_qstates(ek, kmax)
    d = []
    for x = -kmax:kmax
        for y = -kmax:kmax
            for z = -kmax:kmax
                if ek == x * x + y * y + z * z
                    push!(d, Array{Int,1}([x, y, z]))
                end
            end
        end
    end
    d
end

function get_shellindex(ek)
    kmax = ek
    i = sum(collect(get_degen(e, kmax) for e = 1:ek-1))
end


function get_orbindex(k::Array{Int,1})
    println("k² = ", k .* k)
    ek = norm(k .* k)
    i = get_shellindex(ek)
    kstat = get_qstates(ek, ek)
    println("vectors : ", kstat)
    indx = argmax(collect( ([q] ==[k])[1] for q in get_qstates(ek, ek) ))
    println("index : ", indx)
    kvec = kstat[indx]
    # println("kvec  : ", kvec)
    i + indx - 1
end


function main()

    kmax = 10

    # println("Degeneracy of shell 3: ", get_degen(3,kmax))

    # println("Index of shell k*k = 2 : ", get_shellindex(3))

    # println("k-vectors of shell k*k = 2 : ", get_qstates(3,kmax), " => length : ", size(get_qstates(3,kmax),1))

    indx = argmax(collect( ([q] ==[1,1,-1])[1] for q in get_qstates(4, 4) ))

    println("Index of shell k=[-1,1,1] : ")
    println(get_orbindex([-1, 1, 1]))

end

main()
