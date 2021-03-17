using OnlineStats
using DelimitedFiles
using DataFrames
using CSV


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

    e = Ensemble(rs, get_β_internal(θ,N,c), N) # get_β_internal only works for 3D
    updates = Update.([move_particle, add_type_B, remove_type_B, add_type_C, remove_type_C, add_type_D, remove_type_D, add_type_E, remove_type_E, add_remove_kink_chain, shuffle_indices],0,0)#  , add_type_E, remove_type_E, add_remove_kink_chain
                                                                                    #, change_type_B    #


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
    marcov_chain_builders = Array{Task}(undef,Threads.nthreads())# the number of threads set by `julia -t run_threads.jl` (--threads)
    measurements_of_runs = Set{Dict{Symbol,Tuple{OnlineStat,Function}}}()
    for t in 1:Threads.nthreads()
        m = deepcopy(measurements_Mean)
        push!(measurements_of_runs,m)
        marcov_chain_builders[t] = Threads.@spawn(sweep_multithreaded!(NMC, cyc, NEquil, updates, m, e, c))
    end
    for mcb in marcov_chain_builders
        wait(mcb)
    end

    println(" finished.")

    println("parameters:")
    println("θ: ", θ)
    println("rs: ", rs)
    println("N: ", N)

    println("measurements:")
    println("=============")

    # avarage over the uncorrelated mean values of the single runs
    for m in measurements_of_runs
        for (key,(stat,obs)) in m
            if key == :occs
                fit!(first(measurements[key]), mean.(m[:occs][1].stats))
            else
                fit!(first(measurements[key]), mean(stat))
            end
        end
    end

    # print measurements
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

    # print addidtional observables
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
    μT_Ry = E_Ry(μT,λ(e.N,e.rs))
    ΔT_Ry = E_Ry(ΔT,λ(e.N,e.rs))
    μE_Ry = μT_Ry + μWt_Ry
    ΔE_Ry = ΔT_Ry + ΔWt_Ry

    println("W_off_diag", "\t", μW_off_diag, " +/- ", ΔW_off_diag)
    println("W", "\t", μW, " +/- ", ΔW)
    println("E", "\t", μE, " +/- ", ΔE)
    println("W_t_Ry", "\t", μWt_Ry, " +/- ", ΔWt_Ry)
    println("T_Ry", "\t", μT_Ry, " +/- ", ΔT_Ry)

    println("")

    println("occupations:")
    println("============")
    println(mean.(measurements[:occs][1].stats))
    println(std.(measurements[:occs][1].stats))

    # create resultsfile
    # add measurements to file
    df = DataFrame(sign = 1)
    for (k,(f,m)) in measurements
        if typeof(f) == Variance{Float64,Float64,EqualWeight}
            if in(k,[:sign, :K])
                df[!,k] .= mean(f)
                df[!,Symbol(:Δ, k)] .= std(f)/sqrt(Threads.nthreads()-1)
            else
                df[!,k] .= mean(f)/avg_sign
                df[!,Symbol(:Δ, k)] .= std(f)/sqrt(Threads.nthreads()-1)/avg_sign
            end
        end
    end

    # add additional variables to file
    df[!,:W] .= μW
    df[!,:ΔW] .= ΔW
    df[!,:E] .= μE
    df[!,:ΔE] .= ΔE
    df[!,:Wt_Ry] .= μWt_Ry
    df[!,:ΔWt_Ry] .= ΔWt_Ry
    df[!,:T_Ry] .= μT_Ry
    df[!,:ΔT_Ry] .= ΔT_Ry
    df[!,:E_Ry] .= μE_Ry
    df[!,:ΔE_Ry] .= ΔE_Ry

    # create occupation numbers file
    open("test/UEG/Full_CPIMC/out/occNums_$(N)_th$(replace(string(θ),"." => ""))_rs$(replace(string(rs),"." => ""))_steps$((NMC*Threads.nthreads()/cyc)).dat", "w") do io
        writedlm(io, zip(mean.(measurements[:occs][1].stats), std.(measurements[:occs][1].stats)/(NMC*Threads.nthreads()/cyc)))
    end

    CSV.write("test/UEG/Full_CPIMC/out/results_N$(N)_th$(θ)_rs$(rs)_steps$((NMC*Threads.nthreads()/cyc)).csv",df)
end

main()
