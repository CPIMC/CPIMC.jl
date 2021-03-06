using OnlineStats
using DelimitedFiles


include("../../../src/Configuration.jl")
include("../../../src/HEG/model.jl")
include("../../../src/HEG/CPIMC/updates.jl")
include("../../../src/HEG/CPIMC/estimators.jl")
include("../../../src/CPIMC.jl")

function main()
    # MC options
    NMC = 2*10^5
    cyc = 200
    NEquil = 1*10^5

    # system parameters
    θ = 0.125
    rs = 1.75

    S = get_sphere_with_same_spin(OrbitalHEG((0,0,0),1),dk=1)
    #S = union(get_sphere_with_same_spin(OrbitalHEG((0,0,0),1),dk=1),get_sphere_with_same_spin(OrbitalHEG((0,0,0),-1),dk=1))

    #4Particles
    #S = Set{OrbitalHEG{3}}([OrbitalHEG((0,0,0),1), OrbitalHEG((1,0,0),1), OrbitalHEG((0,1,0),1), OrbitalHEG((0,0,1),1)])

    println("#################################################")
    println("N: ", length(S))
    println("θ: ", θ)
    println("rs: ", rs)
    N = length(S)
    c = Configuration(S)

    e = Ensemble(rs, get_β_internal(θ,N,c), N) # get_β_internal only works for 3D

    #We do not need to have change_type_B if we have add_type_C and remove_type_C updates
    updates = [move_particle, add_type_B, remove_type_B ,shuffle_indices, add_type_C, remove_type_C, add_type_D, remove_type_D, add_type_E, remove_type_E, add_remove_kink_chain]

    measurements = Dict(
      :sign => (Variance(), signum)
    , :Ekin => (Variance(), Ekin)
    , :W_diag => (Variance(), W_diag)
    , :K => (Variance(), K)
    , :K_fermion => (Variance(), K)
    , :occs => (Group([Variance() for i in 1:100]), occupations)
    )

    measurements_Mean = Dict(
      :sign => (Mean(), signum)
    , :Ekin => (Mean(), Ekin)
    , :W_diag => (Mean(), W_diag)
    , :K => (Mean(), K)
    , :K_fermion => (Mean(), K)
    , :occs => (Group([Mean() for i in 1:100]), occupations)
    )


    println("Start MC process ... ")
    Marcov_Chain_builders = Array{Task}(undef,Threads.nthreads())#die Anzahl threads ist inital die Anzahl Kerne
    Measurements_of_runs = Set{Dict{Symbol,Tuple{OnlineStat,Function}}}()
    for t in 1:Threads.nthreads()
        m = deepcopy(measurements_Mean)
        push!(Measurements_of_runs,m)
        Marcov_Chain_builders[t] = Threads.@spawn(runMC_multithreaded(NMC, cyc, NEquil, updates, m, e, c))
    end
    for mcb in Marcov_Chain_builders
        wait(mcb)
    end

    println(" finished.")
    println("parameters:")
    println("N: ", length(S))
    println("θ: ", θ)
    println("rs: ", rs)
    println("measurements:")
    println("=============")
    #Avarage over the uncorrolated mean-values of the single runs
    for m in Measurements_of_runs
        for (key,(stat,obs)) in m
            if key == :occs
                fit!(first(measurements[key]), mean.(m[:occs][1].stats))
            else
                fit!(first(measurements[key]), mean(stat))
            end
        end
    end
    #print measurements
    avg_sign = mean(first(measurements[:sign]))
    for (k,(f,m)) in measurements
        if typeof(f) == Variance{Float64,Float64,EqualWeight}
            if in(k,[:sign, :K])
                println(k, "\t", mean(f), " +/- ", std(f)/sqrt(Threads.nthreads()-1))
            else
                println(k, "\t", mean(f)/avg_sign, " +/- ", std(f)/sqrt(Threads.nthreads()-1)/avg_sign)
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
    ΔWt_Ry = E_Ry(ΔW, e::Ensemble)
    μT_Ry = E_Ry(μT,lambda(e.N,e.rs))
    ΔT_Ry = E_Ry(ΔT,lambda(e.N,e.rs))
    #println("W_diag", "\t", μW_diag, " +/- ", ΔW_diag)
    println("W_off_diag", "\t", μW_off_diag, " +/- ", ΔW_off_diag)
    println("W", "\t", μW, " +/- ", ΔW)
    println("E", "\t", μE, " +/- ", ΔE)
    println("W_t_Ry", "\t", μWt_Ry, " +/- ", ΔWt_Ry)
    println("T_Ry", "\t", μT_Ry, " +/- ", ΔT_Ry)


    println("")


    println("occupations:")#plausibilität Prüfen
    println("============")
    println(mean.(measurements[:occs][1].stats))
    println(std.(measurements[:occs][1].stats))

    # create occnumsfile
    #open("test/HEG/rcpimc/out/occNums_$(N)_th$(replace(string(θ),"." => ""))_rs$(replace(string(rs),"." => ""))_Samples$((NMC*Threads.nthreads()/cyc)).dat", "w") do io
    open("CPIMC.jl/test/HEG/rcpimc/out/occNums_$(N)_th$(replace(string(θ),"." => ""))_rs$(replace(string(rs),"." => ""))_Samples$((NMC*Threads.nthreads()/cyc)).dat", "w") do io
           writedlm(io, zip(mean.(measurements[:occs][1].stats), std.(measurements[:occs][1].stats)/(NMC*Threads.nthreads()/cyc)))
    end
    #create resultsfile
    #open("test/HEG/rcpimc/out/results_$(N)_th$(replace(string(θ),"." => ""))_rs$(replace(string(rs),"." => ""))_Samples$((NMC*Threads.nthreads()/cyc)).dat", "w") do io
    open("CPIMC.jl/test/HEG/rcpimc/out/results_$(N)_th$(replace(string(θ),"." => ""))_rs$(replace(string(rs),"." => ""))_Samples$((NMC*Threads.nthreads()/cyc)).dat", "w") do io
        for (k,(f,m)) in measurements
            if typeof(f) == Variance{Float64,Float64,EqualWeight}
                write(io, string(typeof(m).name.mt.name, "\t", mean(f), " +/- ", std(f)/sqrt(Threads.nthreads()-1),"\n"))
                if  (k == :Epot) | (k == :E)
                    write(io, string(typeof(m).name.mt.name,"_t_Ha", "\t", E_Ry(mean(f)-abs_E_mad(e.N, lambda(e.N,e.rs)),lambda(e.N,e.rs))/2, " +/- ", (E_Ry(std(f),lambda(e.N,e.rs))/sqrt(Threads.nthreads()-1)/2),"\n"))
                elseif (k == :Ekin)
                    write(io,string(typeof(m).name.mt.name,"_Ha", "\t", E_Ry(mean(f),lambda(e.N,e.rs))/2, " +/- ", E_Ry(std(f),lambda(e.N,e.rs))/sqrt(Threads.nthreads()-1)/2,"\n"))
                end
            end
        end
    end
end

#Juno.@run(main())
main()
