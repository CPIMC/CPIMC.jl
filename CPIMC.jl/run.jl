using OnlineStats

include("models/canonical/model.jl")
include("models/canonical/updates.jl")
include("MC.jl")
include("System.jl")

function main()
    e = Ensemble(200, 2, 0.1, 50)
    c = Configuration(Set(collect(1:50)))

    updates = Set([move_particle])

    measurements =
    [ (Variance(), totalEnergy)
    , (Group([Variance() for i=1:e.cutoff]), occVec)
    ]

    acc = sweep(10^4, 10, updates, measurements, e, c)

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
#
# updates = Set([move_particle])
# rand(1)
