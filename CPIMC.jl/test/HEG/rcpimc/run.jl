using OnlineStats
using DelimitedFiles


include("../../../src/Configuration.jl")
include("../../../src/HEG/model.jl")
include("../../../src/HEG/RCPIMC/updates.jl")
include("../../../src/HEG/RCPIMC/estimators.jl")
include("../../../src/CPIMC.jl")
include("../../../src/output.jl")


function main()
    # MC options
    NMC = 10^5
    cyc = 20
    NEquil = 10^3
    # system parameters
    θ = 1.0
    rs = 0.5

    S = get_sphere_with_same_spin(OrbitalHEG((0,0,0),0),dk=2) ### use 33 particles
    N = length(S)
    c = Configuration(S)

    println("θ: ", θ)
    println("rs: ", rs)
    println("N: ", N)

    e = Ensemble(rs, get_β_internal(θ,N), N)# get_β_internal only works for 3D
    updates = Update.([move_particle, add_type_B, remove_type_B, change_type_B, shuffle_indices],0,0)

    measurements = Dict(# TODO: type-specification in the construction of the statistic objects (use @code_warntype)
      :Ekin => (Variance(), Ekin)
    , :W_off_diag => (Variance(), W_off_diag)
    , :W_diag => (Variance(), W_diag)
    , :Epot => (Variance(), Epot)
    , :K => (Variance(), K)
    , :E => (Variance(), E)
    , :occs => (Group([Variance(Float64) for i in 1:100]), occupations)# didn't work because of type-specification in Variance(UInt), estimator returned FixedPoint value
    )

    println("Start MC process ... ")
    sweep!(NMC, cyc, NEquil, updates, measurements, e, c)
    println(" finished.")

    print_results(measurements)

    #save_results("out/", measurements, e)
end

#Juno.@run(main())
main()
