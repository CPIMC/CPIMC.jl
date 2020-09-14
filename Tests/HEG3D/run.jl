using OnlineStats

include("../../CPIMC.jl/models/canonical/model.jl")
include("../../CPIMC.jl/models/canonical/updates.jl")
include("../../CPIMC.jl/MC.jl")
include("../../CPIMC.jl/Systems/HEG3D/System.jl")


function main()
    # MC options
    NMC = 10^2
    cyc = 3
    NEquil = 10^1

    # system parameters
    Nb = get_Nb(10)
    N = 4
    theta = 1
    rs = 0.5

    println("Emax = 10, Nb=$Nb, Emax_calculated = $(get_Emax(Nb))")

    e = Ensemble(Nb, rs, get_beta_internal(theta,N), N)
    c = Configuration(Set(collect(1:N)))

    updates = Set([move_particle])

    # for i in 0:10
    #     println("i=$i, Nb=$(get_Nb(i)), Emax=$(get_Emax(get_Nb(i)))")
    # end

    println("cutoff = $(e.cutoff)")

    println("Basis : ", get_basis(e))

    measurements =
    [ (Variance(), Ekin)
    , (Group([Variance() for i=1:e.cutoff]), occVec)
    ]

    # TODO: measurments as Dict()
    # measurements = (:Ekin => (Variance(), Ekin),
    #                 :occN => (Group([Variance() i=1:e.cutoff]), occVec))

    print("Start MC process ... ")

    # sweep(NMC, cyc, NEquil, updates, measurements, e, c)

    println(" finished.")

    println("measurements:")
    println("=============")

    for (f,m) in measurements
        println("measurments $(m) : type=$(typeof(f))")
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
