using DelimitedFiles
using Plots

function EF(N)
    if N == 1
        return 0
    elseif N == 7
        return 1
    elseif N == 19
        return 2
    elseif N == 27
        return 3
    elseif N == 33
        return 4
    else
        println("Warning: Fermi momentum not correct. Please supply information in function EF(N)")
        return -1
    end
end

kF(N) = sqrt(EF(N))

function plot_occNums!(p::Plots.Plot, data, degn, N, Θ, rs; kwargs...)
    x = sqrt.(0:length(data[:,1])-1)# square-root: single-particle energy → absolute of momentum
    y = data[:,1]
    dy = data[:,2]
    ti = .!iszero.(y)# select non-zero results to match with degeneracy_3d.dat
    x = x[ti]
    y = y[ti]
    y = y./degn[1:length(y)]# divide by degeneracy
    dy = dy[ti]
    dy = dy./degn[1:length(y)]# divide by degeneracy

    println("Normalization : sum(n_k) = $(sum(y.*degn[1:length(y)]))")
    println("for $(N) particles.")

    plot(x, y, ribbon=dy, xlabel="k [internal units]", ylabel="n(k)", lw=4, label="n_k")
    vline!([kF(N)], label="kF=$(kF(N))", lw=4, la=0.7)
    title!("Θ=$(Θ), rs=$(rs), N=$(N)")
end
plot_occNums(args...; kwargs...) = plot_occNums!(plot(), args...; kwargs...)

degn = readdlm("data/degeneracy_3d.dat")# load degeneracy

for file in readdir("./out/")
    if occursin("occs_N", file) & occursin(".dat", file)
        n = parse(Int, split(split(file, "occs_N")[2], "_th")[1])
        θ = parse(Float64, split(split(file, "_th")[2], "_rs")[1])
        θ /= 10^floor(log10(θ))
        r = parse(Float64, split(split(file, "_rs")[2], ".dat")[1])
        r /= 10^floor(log10(r))
        data = readdlm(joinpath("out/", file))# load occupation number results
        plot_occNums(data, degn, n, θ, r)
        savefig("out/plots/occs_N$(n)_th$(replace(string(θ), "." => ""))_rs$(replace(string(r), "." => "")).pdf")
    end
end
