using OnlineStats

include("../../src/Configuration.jl")
include("../../src/CPIMC.jl")
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

    S = get_sphere(Orbital((0,0,0),0),dk=2) ### use 19 particles
    N = 33

    println("Number of particles: ", length(S))

    c = Configuration(S)

    e = Ensemble(rs, get_beta_internal(theta,N), N) # get_beta_internal only works for 3D

    updates = [move_particle]

    measurements = Dict(
      :Ekin => (Variance(), Ekin)
    , :occN => (Hist(0:100; left = true, closed = false), occVec)## left border is closed to include the integer values, last right border is open to exclude the right integer
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
    println(measurements[:occN][1])
    println("")
#    println(std.(measurements[:occN][1].stats))## TODO: Variance() of each bin
    println("")
end

main()
