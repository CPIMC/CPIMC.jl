using DelimitedFiles
using Plots

data = readdlm("out/occNums_N33_th10_rs05.dat")

ti = .!iszero.(data[:,2])# select non-zero results to match with degeneracy_3c
x = sqrt.(data[ti,1])#square-root: single-particle energy → absolute of momentum
y = data[ti,2]
y = y./degn[1:length(y)]# divide by degeneracy
y = y./sum(y)# normalize bins

degn = readdlm("data/degeneracy_3d.dat")

plot(x, y, legend=false, xlabel="k", ylabel="n(k)")
