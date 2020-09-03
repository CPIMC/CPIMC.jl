using OnlineStats

include("../../CPIMC.jl/models/canonical/model.jl")
include("../../CPIMC.jl/models/canonical/updates.jl")
include("../../CPIMC.jl/MC.jl")
include("../../CPIMC.jl/Systems/HEG3D/System.jl")


function main()
    # MC options
    NMC = 10^3
    cyc = 4

    # system parameters
    Nb = 200
    N = 50

    e = Ensemble(Nb, 2, 0.1, N)
    c = Configuration(Set(collect(1:N)))

    orblist = get_orblist_UEG(Nb)

    updates = Set([move_particle])

    measurements =
    [ (Variance(), Ekin)
    , (Group([Variance() for i=1:e.cutoff]), occVec)
    ]

    sweep(NMC, cyc, updates, measurements, e, c, orblist)

    println("measurements:")
    println("=============")

    for (f,m) in measurements
        if typeof(f) == Variance{Float64,EqualWeight}
            println(typeof(m).name.mt.name, "\t", mean(f), " +/- ", std(f))
        end
    end

    println("occupations:")
    println("============")
    println(mean.(measurements[2][1].stats))
    println("")
    println(std.(measurements[2][1].stats))
    println("")
end

main()
