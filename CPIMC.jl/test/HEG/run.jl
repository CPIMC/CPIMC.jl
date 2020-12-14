using OnlineStats
using DelimitedFiles


include("../../src/Configuration.jl")
include("../../src/CPIMC.jl")
include("../../src/HEG/model.jl")
include("../../src/HEG/Ideal/updates.jl")
include("../../src/HEG/Ideal/estimators.jl")


function simulation(NMC, cyc, NEquil, theta, rs, dε)
    # MC options
    NMC = NMC
    cyc = cyc
    NEquil = NEquil

    # system parameters
    theta = theta
    rs = rs

    S = get_sphere(Orbital((0,0,0),0),dε=dε)
    N = length(S)

    println("Number of particles: ", N)

    c = Configuration(S)

    e = Ensemble(rs, get_beta_internal(theta,N), N) # get_beta_internal only works for 3D

    updates = [move_particle]

    measurements = Dict(
      :Ekin => (Variance(), Ekin)
    , :occs => (Group([Variance(Int) for i in 1:100]), occupations)
    )

    print("Start MC process ... ")
    runMC(NMC, cyc, NEquil, updates, measurements, e, c)
    println(" finished.")

    print_results(measurements)

    # save_results("./out", measurements, e)
end

function main()
    # MC options
    NMC = 1 * 10^7
    cyc = 10
    NEquil = 10^3
    # system parameters
    T_vals = [1.0]
    R_vals = [0.5]
    dε = 4# 33 particles

    i = 0
    for th in T_vals
        for rs in R_vals
            simulation(NMC, cyc, NEquil, th, rs, dε)
            i += 1
            println(100*i/(length(T_vals)*length(R_vals))," % completed.")
        end
    end
end

main()
