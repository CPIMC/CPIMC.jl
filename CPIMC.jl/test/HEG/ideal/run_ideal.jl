using OnlineStats

include("../../../src/Configuration.jl")
include("../../../src/HEG/model.jl")

include("../../../src/HEG/Ideal/updates.jl")
include("../../../src/HEG/Ideal/estimators.jl")
#include("../../src/HEG/RCPIMC/updates.jl")
#include("../../src/HEG/RCPIMC/estimators.jl")

include("../../../src/CPIMC.jl")

include("../../../src/output.jl")



"""include("CPIMC.jl/src/Configuration.jl")
include("CPIMC.jl/src/CPIMC.jl")
include("CPIMC.jl/src/HEG/model.jl")
include("CPIMC.jl/src/HEG/Ideal/updates.jl")
include("CPIMC.jl/src/HEG/Ideal/estimators.jl")"""

function main()
    # MC options
    NMC = 10^6
    cyc = 10
    NEquil = 10^4

    # system parameters
    θ = 1.0#0.0625
    rs = 0.5#1
    S = get_sphere_with_same_spin(OrbitalHEG((0,0,0),1),dk=2) ### use 33 particles


    println("Number of particles: ", length(S))
    println("θ: ", θ)
    println("rs: ", rs)
    N = length(S)
    println("N: ", N)
    c = Configuration(S)

    e = Ensemble(rs, get_β_internal(θ,N), N) # get_β_internal only works for 3D
    """updates = [move_particle, Add_Type_B, remove_Type_B]

    measurements = Dict(
      :Ekin => (Variance(), Ekin)
    , :W_off_diag => (Variance(), W_off_diag)
    , :W_diag => (Variance(), W_diag)
    , :Epot => (Variance(), Epot)
    , :K => (Variance(), K)
    , :Etot => (Variance(), Etot)
    , :occN => (Group([Variance() for i=1:200]), occVec)
    )"""

    updates = Update.([move_particle],0,0)

    measurements = Dict(
      :Ekin => (Variance(UInt), Ekin)
    , :occs => (Group([Variance(UInt) for i in 1:100]), occupations))

    println("Start MC process ... ")
    sweep!(NMC, cyc, NEquil, updates, measurements, e, c)
    println(" finished.")

    print_results(measurements)

    #save_results("out/", measurements, e)
end

#Juno.@run(main())
main()
