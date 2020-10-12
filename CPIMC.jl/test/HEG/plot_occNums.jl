using DelimitedFiles
using Plots

data = readdlm("out/occNums_N33_th10_rs05.dat")# load occupation number results
degn = readdlm("data/degeneracy_3d.dat")# load degeneracy

ti = .!iszero.(data[:,2])# select non-zero results to match with degeneracy_3c
x = sqrt.(data[ti,1])#square-root: single-particle energy → absolute of momentum
y = data[ti,2]
y = y./degn[1:length(y)]# divide by degeneracy
y = y./sum(y)# normalize bins

plot(x, y, xlabel="k", ylabel="n(k)", lw=4, label="n_k")
vline!([2], label="kF", lw=4, la=0.7)## N=33 ↔ EF=4 → kF=2
title!("Θ=1, rs=0.5, N=33")

savefig("out/plots/occNums_N33_th10_rs05.pdf")
