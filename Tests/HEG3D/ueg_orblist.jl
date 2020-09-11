import DataFrames
import LinearAlgebra: dot

using Plots


"return degeneracy of the shell with energy ek"
function get_degen(ek)
    d = 0
    # TODO: square-root of ek +1 ????
    for x = -ek:ek
        for y = -ek:ek
            for z = -ek:ek
                if ek == x * x + y * y + z * z
                    d += 1
                end
            end
        end
    end
    d
end

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

"return index of shell with energy ek"
function get_shellindex(ek)
    sum(collect(get_degen(e) for e = 0:ek-1))
end

"get index of orbital with momentum vector k"
function get_orbindex(k::Array{Int,1})
    # println("k² = ", k .* k)
    ek = dot(k,k)
    i = get_shellindex(ek)
    # kstat = get_qstates(ek, ek)
    # println("vectors : ", kstat)
    indx = argmax(collect( ([q] == k)[1] for q in get_qstates(ek) ))
    # println("index : ", indx)
    # kvec = kstat[indx]
    # println("kvec  : ", kvec)
    i + indx - 1
end


function main()

    println("--------------------------------------------------------------")
    println("Starting Tests .... ")

    # @time get_shellindex(10)

    k = [rand(-3:3),rand(-3:3),rand(-3:3)]
    ksq = dot(k,k)

    println("Consider k-vector k = $(k), |k|²=$(ksq).")

    si = get_shellindex(ksq)

    println("Shell-Index : ", si)

    println("Degeneracy of this shell : ", get_degen(si))

    println("k-vectors of this shell : ", get_qstates(si))

    println("Orbital Index : ", get_orbindex(k))

    Emax = 10

    # println("Emax = ", Emax)
    print("orblist :\n", DataFrames.DataFrame(indx=0:Emax, qnums=vcat(collect(get_qstates(kk) for kk in 0:Emax)) ) )

end
# sort!()
main()

Emax = 500
plot(1:Emax, collect(@elapsed vcat(collect(get_qstates(kk) for kk in 0:Emax)) for Emax in 1:Emax) , xlabel="Emax", ylabel="Time")
# plot!(x->0.000000002*x^1.5, xlims=(0,500))
