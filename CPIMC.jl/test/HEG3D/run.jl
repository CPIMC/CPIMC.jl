using OnlineStats

include("../../src/CPIMC.jl")
include("../../src/Configuration.jl")
include("../../src/HEG/model.jl")
include("../../src/HEG/Ideal/updates.jl")
include("../../src/HEG/Ideal/estimators.jl")


function main()
    # MC options
    NMC = 5*10^5
    cyc = 3
    NEquil = 10^3

    # system parameters
    theta = 1.0
    rs = 0.5

    S = get_sphere(Orbital((0,0,0),0),dk=2) ### use1 9 particles
    N = 33

    println("Number of particles: ", length(S))

    c = Configuration(S)

    e = Ensemble(rs, get_beta_internal(theta,N), N)

    updates = [move_particle]

    measurements = Dict(
      :Ekin => (Variance(), Ekin)
    , :occN => (Group([Variance() for i=1:200]), occVec)
    )

    print("Start MC process ... ")
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

    println("occupations:")
    println("============")
    println(mean.(measurements[:occN][1].stats))
    println("")
    println(std.(measurements[:occN][1].stats))
    println("")
end

main()
