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

    S = get_sphere(Orbital((0,0,0),0),dk=2) ### use 19 particles
    N = 33

    println("Number of particles: ", length(S))

    c = Configuration(S)

    e = Ensemble(rs, get_beta_internal(theta,N), N) # get_beta_internal only works for 3D

    updates = [move_particle]

    measurements = Dict(
      :Ekin => (Variance(), Ekin)
    , :occs => (Group([Variance() for i in 1:100]), occupations)## left border is closed to include the integer values, last right border is open to exclude the right integer
    , :occ0 => (Variance(), occupation0)
    , :occ0Samples => (CountMap(Int), occupation0)
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
    println(mean.(measurements[:occs][1].stats))
    println(std.(measurements[:occs][1].stats))
    println("occupation of zero-energy state:")
    println("============")
    println(measurements[:occ0Samples])

    # Print to results file
    open("./out/occNums_N$(N)_th$(replace(string(theta),"." => ""))_rs$(replace(string(rs),"." => "")).dat", "w") do io
           writedlm(io, zip(mean.(measurements[:occs][1].stats), std.(measurements[:occs][1].stats)))
       end
end

main()

@time(main())
