using OnlineStats

include("../../src/Configuration.jl")
include("../../src/HEG/model.jl")

include("../../src/HEG/Ideal/updates.jl")
include("../../src/HEG/Ideal/estimators.jl")
#include("../../src/HEG/RCPIMC/updates.jl")
#include("../../src/HEG/RCPIMC/estimators.jl")

include("../../src/CPIMC.jl")


"""include("CPIMC.jl/src/Configuration.jl")
include("CPIMC.jl/src/CPIMC.jl")
include("CPIMC.jl/src/HEG/model.jl")
include("CPIMC.jl/src/HEG/Ideal/updates.jl")
include("CPIMC.jl/src/HEG/Ideal/estimators.jl")"""

function main()
    # MC options
    NMC = 10^6
    cyc = 3
    NEquil = 10^6

    # system parameters
    theta = 0.5
    rs = 0.1

    S = get_orbs_with_spin(get_sphere(Orbital((0,0,0),1),dk=2),1) ### use 33 particles


    println("Number of particles: ", length(S))
    println("theta: ", theta)
    println("rs: ", rs)
    N = length(S)
    println("N: ", N)
    c = Configuration(S)

    e = Ensemble(rs, get_beta_internal(theta,N), N) # get_beta_internal only works for 3D
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

    updates = [move_particle]

    measurements = Dict(
      :Ekin => (Variance(), Ekin)
    )

    println("Start MC process ... ")
    runMC(NMC, cyc, NEquil, updates, measurements, e, c)
    println(" finished.")
    println("measurements:")
    println("=============")

    for (k,(f,m)) in measurements
        if typeof(f) == Variance{Float64,Float64,EqualWeight}
            println(typeof(m).name.mt.name, "\t", mean(f), " +/- ", std(f))
        end
    end

    println("")

    #occupations funktionieren noch nicht fürs WW-System
    """println("occupations:")
    println("============")
    println(mean.(measurements[:occN][1].stats))
    println("")
    println(std.(measurements[:occN][1].stats))
    println("")"""
end

#Juno.@run(main())
main()
