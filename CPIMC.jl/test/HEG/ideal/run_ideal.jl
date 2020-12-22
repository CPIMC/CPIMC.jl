using OnlineStats

include("../../../src/Configuration.jl")
include("../../../src/HEG/model.jl")

include("../../../src/HEG/Ideal/updates.jl")
include("../../../src/HEG/Ideal/estimators.jl")
#include("../../src/HEG/RCPIMC/updates.jl")
#include("../../src/HEG/RCPIMC/estimators.jl")

include("../../../src/CPIMC.jl")


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
    θ = 0.0625
    rs = 1
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

    updates = [move_particle]

    measurements = Dict(
      :Ekin => (Variance(), Ekin)
    , :occs => (Group([Variance() for i in 1:100]), occupations))

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
    println("occupations:")
    println("============")
    println(mean.(measurements[:occs][1].stats))
    println(std.(measurements[:occs][1].stats))

    # Print to results file
    open("CPIMC.jl/test/HEG/rcpimc/out/occNums_N$(N)_th$(replace(string(θ),"." => ""))_rs$(replace(string(rs),"." => "")).dat", "w") do io
           writedlm(io, zip(mean.(measurements[:occs][1].stats), std.(measurements[:occs][1].stats)))
       end
end

#Juno.@run(main())
main()
