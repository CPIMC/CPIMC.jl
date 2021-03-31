using DelimitedFiles
using DataFrames
using CSV

using CPIMC
using CPIMC.Estimators
using CPIMC.UniformElectronGas

import CPIMC: 

function measure(measurements, e, c)
    s = signum(c)
    for (key,(stat,obs)) in measurements
        if in(key,[:sign, :K, :T1c])
            fit!(stat, obs(e,c))
        else
            fit!(stat, obs(e,c)*s)
        end
    end
end

function W_off_diag(e::Ensemble, avg_K::Float64)
    return (-(avg_K/e.β))
end

function abs_E_madelung(N::Int, λ::Float64) #internal units
    return 2.83729747948527 * pi/2.0 * N * λ
end

function E_int_from_Hartree(E_Ha::Float64, λ::Float64)
    return (E_Ha /(16/((2*pi)^4 * (λ/2)^2) * 0.5))
end

function E_int_from_Hartree(E_Ha::Float64, e::Ensemble)
    return (E_Ha /(16/((2*pi)^4 * (e.λ/2)^2) * 0.5))
end

function E_Ry(E_internal::Float64, λ::Float64)
    return (E_internal * 16/((2*pi)^4 * λ^2))
end

function E_Ry(E_internal::Float64, e::Ensemble)
    return (E_internal * 16/((2*pi)^4 * e.λ^2))
end

function E_Ha(E_internal::Float64, λ::Float64)
    return (E_internal * 16/((2*pi)^4 * λ^2) * 0.5)
end

function E_Ha(E_internal::Float64, e::Ensemble)
    return (E_internal * 16/((2*pi)^4 * e.λ^2) * 0.5)
end


function Et_Ry(E_internal::Float64, e::Ensemble)
    return (E_Ry(E_internal-abs_E_madelung(e.N, e.λ),e.λ))
end

function Et_Ha(E_internal::Float64, e::Ensemble)
    return E_Ha(E_internal-abs_E_madelung(e.N, e.λ),e.λ)
end


function sweep_multithreaded!(steps::Int, sampleEvery::Int, throwAway::Int, updates::Array{Update,1}, measurements, e::Ensemble, c::Configuration; kwargs...)
    @assert(length(c.kinks) == 0)
    c = Configuration(copy(c.occupations))# c should be a different object for each thread
    " equilibration "
    if (Threads.threadid() == 1)
        println("\nstarting equilibration")
    end
    k = 1# progress counter
    for i in 1:throwAway
        if (i%(throwAway/100) == 0)
            println("                 "^(Threads.threadid()-1),"T",Threads.threadid(), " eq: ",k,"/100","; K: ",length(c.kinks))

            k+=1
        end
        #TODO Use reentrantlook for Update counters?
        update!(c, e, updates; kwargs...)
    end
    if (Threads.threadid() == 1)
        println("\nstarting Simulation")
    end
    i = 0
    k = 1#print progress
    #global add_E_counter = 0
    #global remove_E_counter = 0
    while i < steps
        #print progress
        if (i%(steps/100) == 0) #& (Threads.threadid() == 1)
            println("                 "^(Threads.threadid()-1),"T",Threads.threadid(), " ",k,"/100","; K: ",length(c.kinks))
            k+=1
        end

        " MC step "
        update!(c, e, updates; kwargs...)

        "measurement"
        if i % sampleEvery == 0
            " calculate observables "
            measure(measurements, e, c)
        end
        i += 1
    end
    println("\nThread",Threads.threadid(),"finished")
end


const ex_radius = 3 # maximum radius for exitation


#To run on a Linux System use "julia --threads NT run_threads.jl", where NT is the
#desired number of Threads.
#Inside the Code you can use Threads.nthreads() to check how many Threads
#the Programm uses.
function main()
    # MC options
    NMC = 5*10^5
    cyc = 50
    N_Runs = 24
    NEquil = 10^5
    # system parameters
    θ = 0.125
    rs = 2.0

    # use 7 particles
    S = sphere_with_same_spin(OrbitalHEG((0,0,0),Up),dk=1)
    #S = sphere(OrbitalHEG((0,0,0),Up),dk=2)
    N = length(S)
    c = Configuration(S)

    println("parameters:")
    println("θ: ", θ)
    println("rs: ", rs)
    println("N: ", N)

    e = CEnsemble(λ(N,rs,d), β(θ,N,fractional_spin_polarization(c)), N) # get_β_internal only works for 3D
    updates = Update.([move_particle, add_type_B, remove_type_B, add_type_C, remove_type_C, add_type_D, remove_type_D, add_type_E, remove_type_E, add_remove_kink_chain, shuffle_indices],0,0,0)#, change_type_B


    measurements = Dict(
      :sign => (Variance(), signum)
    , :Ekin => (Variance(), Ekin)
    , :W_diag => (Variance(), W_diag)
    , :K => (Variance(), K)
    , :K_fermion => (Variance(), K)
    , :T1c => (Variance(), right_type_1_count)
    , :occs => (Group([Variance() for i in 1:100]), occupations)
    )

    measurements_Mean = Dict(
      :sign => (Mean(), signum)
    , :Ekin => (Mean(), Ekin)
    , :W_diag => (Mean(), W_diag)
    , :K => (Mean(), K)
    , :K_fermion => (Mean(), K)
    , :T1c => (Mean(), right_type_1_count)
    , :occs => (Group([Mean() for i in 1:100]), occupations)
    )

    println("Start MC process ... ")
    measurements_of_runs = Set{Dict{Symbol,Tuple{OnlineStat,Function}}}()


    Threads.@threads for t in 1:N_Runs
        m = deepcopy(measurements_Mean)
        push!(measurements_of_runs,m)
        sweep_multithreaded!(NMC, cyc, NEquil, updates, m, e, c)
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
            if in(k,[:sign, :K, :T1c])
                println(k, "\t", mean(f), " +/- ", std(f)/sqrt(N_Runs-1))
            else
                println(k, "\t", mean(f)/avg_sign, " +/- ", std(f)/sqrt(N_Runs-1)/avg_sign)
            end
        end
    end

    # print addidtional observables
    μW_diag = mean(first(measurements[:W_diag]))/avg_sign
    ΔW_diag = std(first(measurements[:W_diag]))/sqrt(N_Runs-1)/avg_sign
    μW_off_diag = W_off_diag(e::Ensemble, mean(first(measurements[:K_fermion]))/avg_sign)
    ΔW_off_diag = abs(W_off_diag(e::Ensemble, std(first(measurements[:K_fermion]))/sqrt(N_Runs-1)/avg_sign))
    μT = mean(first(measurements[:Ekin]))/avg_sign
    ΔT = std(first(measurements[:Ekin]))/sqrt(N_Runs-1)/avg_sign
    μW = μW_diag + μW_off_diag
    ΔW = ΔW_diag + ΔW_off_diag
    μE = μW + μT
    ΔE = ΔW + ΔT
    μWt_Ha = Et_Ha(μW, e::Ensemble)
    ΔWt_Ha = E_Ha(ΔW, e::Ensemble)
    μT_Ha = E_Ha(μT,e.λ)
    ΔT_Ha = E_Ha(ΔT,.λ)
    μE_Ha = μT_Ha + μWt_Ha
    ΔE_Ha = ΔT_Ha + ΔWt_Ha

    println("W_off_diag", "\t", μW_off_diag, " +/- ", ΔW_off_diag)
    println("W", "\t", μW, " +/- ", ΔW)
    println("E", "\t", μE, " +/- ", ΔE)
    println("W_t_Ha", "\t", μWt_Ha, " +/- ", ΔWt_Ha)
    println("T_Ha", "\t", μT_Ha, " +/- ", ΔT_Ha)

    println("")

    println("occupations:")
    println("============")
    println(mean.(measurements[:occs][1].stats))
    println(std.(measurements[:occs][1].stats))

    println("acceptance ratios:")
    println("============")
    for u in updates
        println("$(u.update):\t$(u.proposed) proposed,\t$(u.accepted) accepted,\t$(u.trivial) trivial,\tratio(acc/prop) : $(u.accepted/u.proposed), ratio(acc/(prop-triv)) : $(u.accepted/(u.proposed-u.trivial))")
    end


    # create resultsfile
    # add measurements to file
    df = DataFrame(sign = 1)
    for (k,(f,m)) in measurements
        if typeof(f) == Variance{Float64,Float64,EqualWeight}
            if in(k,[:sign, :K])
                df[!,k] .= mean(f)
                df[!,Symbol(:Δ, k)] .= std(f)/sqrt(N_Runs-1)
            else
                df[!,k] .= mean(f)/avg_sign
                df[!,Symbol(:Δ, k)] .= std(f)/sqrt(N_Runs-1)/avg_sign
            end
        end
    end

    # add additional variables to file
    df[!,:W] .= μW
    df[!,:ΔW] .= ΔW
    df[!,:E] .= μE
    df[!,:ΔE] .= ΔE
    df[!,:Wt_Ha] .= μWt_Ha
    df[!,:ΔWt_Ha] .= ΔWt_Ha
    df[!,:T_Ha] .= μT_Ha
    df[!,:ΔT_Ha] .= ΔT_Ha
    df[!,:E_Ha] .= μE_Ha
    df[!,:ΔE_Ha] .= ΔE_Ha

    # create occupation numbers file
    open("test/UEG/Full_CPIMC/out/occNums_$(N)_th$(replace(string(θ),"." => ""))_rs$(replace(string(rs),"." => ""))_steps$((NMC*N_Runs/cyc)).dat", "w") do io
        writedlm(io, zip(mean.(measurements[:occs][1].stats), std.(measurements[:occs][1].stats)/(NMC*N_Runs/cyc)))
    end

    CSV.write("test/UEG/Full_CPIMC/out/results_N$(N)_th$(θ)_rs$(rs)_steps$((NMC*N_Runs/cyc)).csv",df)
end

main()

