#] activate .
using DelimitedFiles
using DataFrames
using CSV
using Revise
using CPIMC
using CPIMC.Estimators
using CPIMC.PlaneWaves
using CPIMC.UniformElectronGas
using StaticArrays


using OnlineStats

import CPIMC: move_particle, add_type_B, remove_type_B, change_type_B, add_type_C, remove_type_C, add_type_D, remove_type_D, add_type_E, remove_type_E, add_remove_kink_chain, shuffle_indices, equlibrate_diagonal!

import CPIMC: measure!, update!

W_off_diag(e::Ensemble, avg_K::Float64) = -avg_K/e.β
abs_E_madelung(N::Int, λ::Float64) = 2.83729747948527 * pi/2.0 * N * λ
E_int_from_Hartree(E_Ha::Float64, λ::Float64) = E_Ha /(16/((2*pi)^4 * (λ/2)^2) * 0.5)
E_int_from_Hartree(E_Ha::Float64, e::Ensemble) = E_Ha /(16/((2*pi)^4 * (e.λ/2)^2) * 0.5)
E_Ry(E_internal::Float64, λ::Float64) = E_internal * 16/((2*pi)^4 * λ^2)
E_Ry(E_internal::Float64, e::Ensemble) = E_internal * 16/((2*pi)^4 * e.λ^2)
E_Ha(E_internal::Float64, λ::Float64) = E_internal * 16/((2*pi)^4 * λ^2) * 0.5
E_Ha(E_internal::Float64, e::Ensemble) = E_internal * 16/((2*pi)^4 * e.λ^2) * 0.5
Et_Ry(E_internal::Float64, e::Ensemble) = E_Ry(E_internal-abs_E_madelung(e.N, e.λ),e.λ)
Et_Ha(E_internal::Float64, e::Ensemble) = E_Ha(E_internal-abs_E_madelung(e.N, e.λ),e.λ)

#To run on a Linux System use "julia --threads NT run_threads.jl", where NT is the
#desired number of Threads.
#Inside the Code you can use Threads.nthreads() to check how many Threads
#the Programm uses.
function main()
    # MC options
    NMC = 10^5
    cyc = 50
    N_Runs = 24
    NEquil = 2*10^4
    # system parameters
    θ = 0.125
    rs = 2.0

    # use 7 particles
    S = sphere_with_same_spin(PlaneWave((0,0,0),Up),dk=1)
    #S = sphere(PlaneWave((0,0,0),Up),dk=2)
    N = length(S)
    c = Configuration(S)
    d = dimension(c.occupations)
    ξ = fractional_spin_polarization(S)

    println("parameters:")
    println("θ: ", θ)
    println("rs: ", rs)
    println("N: ", N)
    println("ξ: ", ξ)

    e = CEnsemble(λ(N,rs,d), β(θ,N,ξ,d), N)


    update_names = [move_particle, add_type_B, remove_type_B, add_type_C, remove_type_C, add_type_D, remove_type_D, add_type_E, remove_type_E, add_remove_kink_chain, shuffle_indices]
    update_counters = []
    for _ in 1:length(update_names)
        push!(update_counters, UpdateCounter())
    end


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
        me = deepcopy(measurements_Mean)
        updates = Array{Tuple{Function, UpdateCounter},1}()
        for up_name in update_names
            push!(updates, (up_name,UpdateCounter()))
        end
        push!(measurements_of_runs,me)
        c_new = Configuration(copy(c.occupations))
            equlibrate_diagonal!(UEG(), e, c_new)
        sweep!(UEG(), e, c_new, updates, me, NMC, cyc, NEquil)
        lock(ReentrantLock()) do
            for i in 1:length(updates)
                update_counters[i] += updates[i][2]
            end
        end
    end
    #Addup counters


    println(" finished.")

    println("parameters:")
    println("θ: ", θ)
    println("rs: ", rs)
    println("N: ", N)

    println("measurements:")
    println("=============")

    # average over the uncorrelated mean values of the single runs
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
    μWt_Ha = Et_Ha(μW, e)
    ΔWt_Ha = E_Ha(ΔW, e)
    μT_Ha = E_Ha(μT,e.λ)
    ΔT_Ha = E_Ha(ΔT,e.λ)
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
    for i in 1:length(update_counters)
        uc = update_counters[i]
        up = update_names[i]
        println("$(up):\t$(uc.proposed) proposed,\t$(uc.accepted) accepted,\t$(uc.trivial) trivial,\tratio(acc/prop) : $(uc.accepted/uc.proposed), ratio(acc/(prop-triv)) : $(uc.accepted/(uc.proposed-uc.trivial))")
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
    open("examples/out/occNums_$(N)_th$(replace(string(θ),"." => ""))_rs$(replace(string(rs),"." => ""))_steps$((NMC*N_Runs/cyc)).dat", "w") do io
        writedlm(io, zip(mean.(measurements[:occs][1].stats), std.(measurements[:occs][1].stats)/(NMC*N_Runs/cyc)))
    end

    CSV.write("examples/out/results_N$(N)_th$(θ)_rs$(rs)_steps$((NMC*N_Runs/cyc)).csv",df)
end

main()
