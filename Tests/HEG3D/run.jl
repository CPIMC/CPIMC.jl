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
    Ec = 10
    Nb = get_Nb(Ec)
    N = 4
    theta = 1
    rs = 0.5

    e = Ensemble(Nb, rs, get_beta_internal(theta,N), N)

    # "return a Set{Orbital} of all Orbitals in the Basis"
    # function get_basis(Emax)
    # end
    # TODO: use empty constructor
    # o = Set{Orbital{3}}([Orbital{3}(SVector{3, Int16}(0,0,0))])
    # o = Orbital{3}(SVector{3,Int16}(1,2,3))
    o = Set{Orbital{3}}()
    # Orbital{3}( SVector{3, Int16}(convert(Int16,1),convert(Int16,2),convert(Int16,3)) )
    # kc = Int(floor(sqrt(Ec)))
    #

    # print(o)
    # c = Configuration(Set(Orbital.(SVector{3,Int16}.(rand(Int16,10),rand(Int16,10),rand(Int16,10)))))
    # c = Configuration(Set(collect(1:N)))
    # b = Basis(o)
    # println("random orbital : ", rand(b.Orbitals) )
    # println("2 random orbitals : ", [rand(b.Orbitals) for i in 1:4] )
    # c = Configuration(ran)
    # updates = Set([move_particle])
    # for i in 0:10
    #     println("i=$i, Nb=$(get_Nb(i)), Emax=$(get_Emax(get_Nb(i)))")
    # end
    # println("Basis : ", get_basis(e))

    measurements =
    [ (Variance(), Ekin)
    , (Group([Variance() for i=1:e.cutoff]), occVec)
    ]
    # TODO: measurments as Dict()
    # measurements = (:Ekin => (Variance(), Ekin),
    #                 :occN => (Group([Variance() i=1:e.cutoff]), occVec))
    print("Start MC process ... ")
    # sweep(NMC, cyc, NEquil, updates, measurements, e, c, o)
    println(" finished.")
    # println("measurements:")
    # println("=============")
    #
    # for (f,m) in measurements
    #     println("measurments $(m) : type=$(typeof(f))")
    #     if typeof(f) == Variance{Float64,EqualWeight}
    #         println(typeof(m).name.mt.name, "\t", mean(f), " +/- ", std(f))
    #     end
    # end
    #
    # println("occupations:")
    # println("============")
    # println(mean.(measurements[2][1].stats))
    # println("")
    # println(std.(measurements[2][1].stats))
    # println("")
end

main()
