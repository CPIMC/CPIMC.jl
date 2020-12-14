using OnlineStats
using DelimitedFiles

include("../../src/HEG/model.jl")
include("../../src/Configuration.jl")
include("../../src/HEG/Ideal/estimators.jl")
include("../../src/CPIMC.jl")
include("../../src/HEG/Ideal/updates.jl")
include("../../src/output.jl")


function simulation(NMC, cyc, NEquil, theta, rs, dε)

    S = get_sphere(Orbital((0,0,0),0),dε=dε)
    N = length(S)
    println("Number of particles: ", N)

    c = Configuration(S)

    e = Ensemble(rs, get_beta_internal(theta,N), N) # get_beta_internal only works for 3D

    updates = Update.([move_particle!], 0, 0)

    measurements = Dict(
        :Ekin => (Variance(), Ekin)
    ,   :occs => (Group([Variance(UInt) for i = 1:100]), occupations)
    )

    print("Start MC process ... ")
    sweep(NMC, cyc, NEquil, updates, measurements, e, c)
    println(" finished.")

    print_results(measurements)

    save_results("/home/vorlautv/git/cpimc2020/CPIMC.jl/test/HEG/out", measurements, e)
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

"perform an MC run without taking measurements, thus profiling the MC process only"
function profile_MC()
    # MC options: perform only equlibration
    NMC = 0
    cyc = 1
    NEquil = 10^7
    # system parameters
    theta = 1.0
    rs = 0.5
    dε = 4# 33 particles
    S = get_sphere(Orbital((0,0,0),0),dε=dε)
    N = length(S)
    c = Configuration(S)
    e = Ensemble(rs, get_beta_internal(theta,N), N) # get_beta_internal only works for 3D
    updates = Update.([move_particle!], 0, 0)

    measurements = Dict(
        :Ekin => (Variance(), Ekin)
    ,   :occs => (Group([Variance(UInt) for i = 1:100]), occupations)
    )

    sweep(NMC, cyc, NEquil, updates, measurements, e, c)
end

main()
