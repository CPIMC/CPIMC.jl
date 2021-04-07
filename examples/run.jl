#] activate .
using StaticArrays
using OnlineStats
using DelimitedFiles
using Revise
using CPIMC
using CPIMC.PlaneWaves
using CPIMC.Estimators
using CPIMC.UniformElectronGas
import CPIMC: move_particle, add_type_B, remove_type_B, add_type_C, remove_type_C, add_type_D, remove_type_D, add_type_E, remove_type_E, add_remove_kink_chain, shuffle_indices

function main()
    # MC options
    NMC = 3 * 10^5
    cyc = 50
    NEquil = 5*10^4
    # system parameters
    θ = 0.125
    rs = 2.0

    # use 7 particles
    S = sphere_with_same_spin(PlaneWave((0,0,0)),dk=1)
    N = length(S)
    ξ = fractional_spin_polarization(S)
    c = Configuration(S)
    d = dimension(c.occupations)

    println("parameters:")
    println("θ: ", θ)
    println("rs: ", rs)
    println("N: ", N)
    println("ξ: ", ξ)

    e = CEnsemble(λ(N, rs, d), β(θ, N, ξ, d), N)


    update_names = [move_particle, add_type_B, remove_type_B, add_type_C, remove_type_C, add_type_D, remove_type_D, add_type_E, remove_type_E, add_remove_kink_chain, shuffle_indices]
    updates = Array{Tuple{Function,MArray{Tuple{3},Int64,1,3}},1}()
    for up_name in update_names
        push!(updates, (up_name,MVector((0,0,0))))
    end

    measurements = Dict(# TODO: type-specification in the construction of the statistic objects (use @code_warntype)
      :sign => (Variance(), signum)
    , :Ekin => (Variance(), Ekin)
    , :W_diag => (Variance(), W_diag)
    , :K => (Variance(), K)
    , :K_fermion => (Variance(), K)
    , :occs => (Group([Variance(Float64) for i in 1:100]), occupations)
    )

    println("Start MC process ... ")
    sweep!(UEG(), e, c, updates, measurements, NMC, cyc, NEquil)
    println(" finished.")

    println("parameters:")
    println("θ: ", θ)
    println("rs: ", rs)
    println("N: ", N)

    for u in updates
        println("$(u[1]):\t$(u[2][1]) proposed,\t$(u[2][2]) accepted,\t$(u[2][3]) trivial,\tratio(acc/prop) : $(u[2][2]/u[2][1]), ratio(acc/(prop-triv)) : $(u[2][2]/(u[2][1]-u[2][3]))")
    end


    print_results(measurements, e)
end

main()
#Juno.@run(main())
