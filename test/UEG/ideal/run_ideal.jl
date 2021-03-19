using OnlineStats

include("../../../src/Configuration.jl")
include("../../../src/UEG/model.jl")
include("../../../src/Updates/Ideal-Updates.jl")
include("../../../src/UEG/estimators.jl")
include("../../../src/CPIMC.jl")
include("../../../src/output.jl")


function main()
    # MC options
    NMC = 10^6
    cyc = 10
    NEquil = 10^4
    # system parameters
    θ = 1.0
    rs = 0.5

    S = sphere_with_same_spin(OrbitalHEG((0,0,0)),dk=2)### use 33 particles
    N = length(S)
    ξ = fractional_spin_polarization(S)
    c = Configuration(S)

    println("θ: ", θ)
    println("rs: ", rs)
    println("N: ", N)

    e = Ensemble(rs, β(θ,N,ξ), N) # β only works for 3D

    updates = Update.([move_particle],0,0,0)

    measurements = Dict(
      :Ekin => (Variance(UInt), Ekin)
    , :occs => (Group([Variance(UInt) for i in 1:100]), occupations))

    println("Start MC process ... ")
    sweep!(NMC, cyc, NEquil, updates, measurements, e, c)
    println(" finished.")

    print_results(measurements, e)

    #save_results("out/", measurements, e)
end

#Juno.@run(main())
main()
