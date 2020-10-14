using DelimitedFiles
using Plots

data = readdlm("out/occNums_N33_th10_rs05.dat")# load occupation number results
degn = readdlm("data/degeneracy_3d.dat")# load degeneracy

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

plot(x, y, ribbon=dy, xlabel="k", ylabel="n(k)", lw=4, label="n_k", ylims=(0,1))
vline!([2], label="kF", lw=4, la=0.7)## N=33 ↔ EF=4 → kF=2
title!("Θ=1, rs=0.5, N=33")

savefig("out/plots/occNums_N33_th10_rs05.pdf")
