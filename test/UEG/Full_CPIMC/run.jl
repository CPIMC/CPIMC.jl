using OnlineStats
using DelimitedFiles


include("../../../src/Configuration.jl")
include("../../../src/UEG/model.jl")
include("../../../src/Updates/Ideal-Updates.jl")
include("../../../src/Updates/Other-Updates.jl")
include("../../../src/Updates/Type-A-Updates.jl")
include("../../../src/Updates/Type-B-Updates.jl")
include("../../../src/Updates/Type-C-Updates.jl")
include("../../../src/Updates/Type-D-Updates.jl")
include("../../../src/Updates/Type-E-Updates.jl")
include("../../../src/UEG/estimators.jl")
include("../../../src/CPIMC.jl")
include("../../../src/output.jl")


const ex_radius = 3 # maximum radius for exitation


function main()
    # MC options
    NMC = 5 * 10^5
    cyc = 50
    NEquil = 10^5
    # system parameters
    θ = 0.125
    rs = 2.0

    # use 7 particles
    S = get_sphere_with_same_spin(OrbitalHEG((0,0,0),1),dk=1)
    N = length(S)
    c = Configuration(S)

    println("parameters:")
    println("θ: ", θ)
    println("rs: ", rs)
    println("N: ", N)

    e = Ensemble(rs, get_β_internal(θ,N,c), N) # TODO: get_β_internal only works for 3D
    updates = Update.([move_particle, add_type_B, remove_type_B, add_type_C, remove_type_C, add_type_D, remove_type_D, add_type_E, remove_type_E, add_remove_kink_chain, shuffle_indices],0,0) # change_type_B is not required for CPIMC, but for RCPIMC

    measurements = Dict(# TODO: type-specification in the construction of the statistic objects (use @code_warntype)
      :sign => (Variance(), signum)
    , :Ekin => (Variance(), Ekin)
    , :W_diag => (Variance(), W_diag)
    , :K => (Variance(), K)
    , :K_fermion => (Variance(), K)
    , :occs => (Group([Variance(Float64) for i in 1:100]), occupations)
    )

    println("Start MC process ... ")
    sweep!(NMC, cyc, NEquil, updates, measurements, e, c)
    println(" finished.")

    println("parameters:")
    println("θ: ", θ)
    println("rs: ", rs)
    println("N: ", N)

    print_results(measurements)
    #save_results("out/", measurements, e)
end

#Juno.@run(main())
main()
