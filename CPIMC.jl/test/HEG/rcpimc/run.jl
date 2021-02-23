using OnlineStats
using DelimitedFiles


include("../../../src/Configuration.jl")
include("../../../src/HEG/model.jl")
include("../../../src/HEG/RCPIMC/updates.jl")
include("../../../src/HEG/RCPIMC/estimators.jl")
include("../../../src/CPIMC.jl")


"""include("CPIMC.jl/src/Configuration.jl")
include("CPIMC.jl/src/HEG/model.jl")
include("CPIMC.jl/src/CPIMC.jl")
include("CPIMC.jl/src/HEG/RCPIMC/updates.jl")
include("CPIMC.jl/src/HEG/RCPIMC/estimators.jl")"""

function main()
    # MC options
    NMC = 10^5###############################################################
    cyc = 200
    NEquil = 10*10^5
    #auffälligerBalken um schwer übersehbaren unterschied im vergleich zu run_threads herzustellen
    """#####################################################################
    ########################################################################
    ####################################################################"""
    # system parameters
    θ = 0.125
    rs = 1.5

    S = get_sphere_with_same_spin(OrbitalHEG((0,0,0),1),dk=1)

    #2Particles
    #S = Set{OrbitalHEG{3}}([OrbitalHEG((0,0,0),1), OrbitalHEG((1,0,0),1), OrbitalHEG((0,1,0),1), OrbitalHEG((0,0,1),1)])

    println("Number of particles: ", length(S))
    println("θ: ", θ)
    println("rs: ", rs)
    N = length(S)
    println("N: ", N)
    c = Configuration(S)

    e = Ensemble(rs, get_β_internal(θ,N), N) # get_β_internal only works for 3D
    updates = [move_particle, add_type_B, remove_type_B ,shuffle_indices, add_type_C, remove_type_C, add_type_D, remove_type_D, add_type_E, remove_type_E, Add_remove_Kink_Chain] #, change_type_B

    measurements = Dict(
      :sign => (Variance(), signum)
    , :Ekin => (Variance(), Ekin)
    , :W_diag => (Variance(), W_diag)
    , :K => (Variance(), K)
    , :K_fermion => (Variance(), K)
    , :occs => (Group([Variance() for i in 1:100]), occupations)
    )

    println("Start MC process ... ")
    runMC(NMC, cyc, NEquil, updates, measurements, e, c)
    println(" finished.")
    println("parameters:")
    println("N: ", length(S))
    println("θ: ", θ)
    println("rs: ", rs)
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

    #occupations funktionieren noch nicht fürs WW-System
    println("occupations:")
    println("============")
    println(mean.(measurements[:occs][1].stats))
    println(std.(measurements[:occs][1].stats))

    # Print to results file
    #open("../out/occNums_N$(N)_th$(replace(string(θ),"." => ""))_rs$(replace(string(rs),"." => "")).dat", "w") do io
    #       writedlm(io, zip(mean.(measurements[:occs][1].stats), std.(measurements[:occs][1].stats)))
    #   end
end

#Juno.@run(main())
main()
