using OnlineStats
using DelimitedFiles


include("../../src/Configuration.jl")
include("../../src/CPIMC.jl")
include("../../src/HEG/model.jl")
include("../../src/HEG/Ideal/updates.jl")
include("../../src/HEG/Ideal/estimators.jl")


function main()
    # MC options
    NMC = 1*10^6
    cyc = 10
    NEquil = 10^3

    # system parameters
    theta = 1.0
    rs = 0.5

    S = get_sphere(Orbital((0,0,0),0),dk=2)
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

    save_results("./out", measurements, e)
end

main()
