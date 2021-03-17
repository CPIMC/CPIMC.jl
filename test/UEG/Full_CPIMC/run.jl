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
const ex_radius = 3 #max Radius for exitation
function main()
    # MC options
    NMC = 5*10^6
    cyc = 50
    NEquil = 5*10^5
    # system parameters
    θ = 0.5
    rs = 0.5

    #S = union!(get_sphere_with_same_spin(OrbitalHEG((0,0,0),1),dk=1), get_sphere_with_same_spin(OrbitalHEG((0,0,0),-1),dk=1))
    S = get_sphere_with_same_spin(OrbitalHEG((0,0,0),1),dk=2)
    N = length(S)
    c = Configuration(S)

    println("θ: ", θ)
    println("rs: ", rs)
    println("N: ", N)

    # e = Ensemble(rs, get_β_internal(θ,N), N)# get_β_internal only works for 3D
    # updates = Update.([move_particle, add_type_B, remove_type_B, change_type_B, shuffle_indices],0,0)

    e = Ensemble(rs, get_β_internal(θ,N,c), N) # get_β_internal only works for 3D
    updates = Update.([move_particle, add_type_B, remove_type_B, change_type_B, shuffle_indices, add_type_C, remove_type_C, add_type_D, remove_type_D, add_type_E, remove_type_E, add_remove_kink_chain],0,0)#, change_type_B  , add_type_C, remove_type_C, add_type_D, remove_type_D, add_type_E, remove_type_E, add_remove_kink_chain

    measurements = Dict(# TODO: type-specification in the construction of the statistic objects (use @code_warntype)
      :sign => (Variance(), signum)
    , :Ekin => (Variance(), Ekin)
    , :W_diag => (Variance(), W_diag)
    , :K_fermion => (Variance(), K)
    , :K => (Variance(), K)
    , :occs => (Group([Variance(Float64) for i in 1:100]), occupations)# didn't work because of type-specification in Variance(UInt), estimator returned FixedPoint value
    )

    println("Start MC process ... ")
    sweep!(NMC, cyc, NEquil, updates, measurements, e, c)
    println(" finished.")
    # println("parameters:")
    # println("N: ", length(S))
    # println("θ: ", θ)
    # println("rs: ", rs)
    println("measurements:")
    println("=============")

    avg_sign = mean(first(measurements[:sign]))
    for (k,(f,m)) in measurements
        if typeof(f) == Variance{Float64,Float64,EqualWeight}
            if in(k,[:sign, :K])
                println(k, "\t", mean(f), " +/- ", std(f)/sqrt(NMC/cyc))
            else
                println(k, "\t", mean(f)/avg_sign, " +/- ", std(f)/sqrt(NMC/cyc)/avg_sign)
            end
        end
    end

    #print addidtional observables
    μW_diag = mean(first(measurements[:W_diag]))/avg_sign
    ΔW_diag = std(first(measurements[:W_diag]))/sqrt(Threads.nthreads()-1)/avg_sign
    μW_off_diag = W_off_diag(e::Ensemble, mean(first(measurements[:K_fermion]))/avg_sign)
    ΔW_off_diag = abs(W_off_diag(e::Ensemble, std(first(measurements[:K_fermion]))/sqrt(Threads.nthreads()-1)/avg_sign))
    μT = mean(first(measurements[:Ekin]))/avg_sign
    ΔT = std(first(measurements[:Ekin]))/sqrt(Threads.nthreads()-1)/avg_sign
    μW = μW_diag + μW_off_diag
    ΔW = ΔW_diag + ΔW_off_diag
    μE = μW + μT
    ΔE = ΔW + ΔT
    μWt_Ry = Et_Ry(μW, e::Ensemble)
    ΔWt_Ry = Et_Ry(ΔW, e::Ensemble)
    μT_Ry = E_Ry(μT,lambda(e.N,e.rs))
    ΔT_Ry = E_Ry(ΔT,lambda(e.N,e.rs))
    println("W_diag", "\t", μW_diag, " +/- ", ΔW_diag)
    println("W_off_diag", "\t", μW_off_diag, " +/- ", ΔW_off_diag)
    println("W", "\t", μW, " +/- ", ΔW)
    println("E", "\t", μE, " +/- ", ΔE)
    println("W_t_Ry", "\t", Et_Ry(μW,e), " +/- ", Et_Ry(ΔW,e))
    println("T_Ry", "\t", E_Ry(μT,e), " +/- ", E_Ry(ΔT,e))
    println("")

    #print_results(measurements)
    #save_results("out/", measurements, e)
end

#Juno.@run(main())
@time main()
