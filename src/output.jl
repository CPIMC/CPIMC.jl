using DelimitedFiles

" print all measured observables and calculate further quantities "
function print_results(measurements, e::Ensemble)

    println("measurements:")
    println("=============")

    # calculate average sign
    if haskey(measurements, :sign)
        avg_sign = mean(first(measurements[:sign]))
    else
        avg_sign = 1.0
    end

    for (k,(f,m)) in measurements
        if typeof(f) == Variance{Float64,Float64,EqualWeight}
            if in(k,[:sign, :K])
                println(k, "\t", mean(f), " +/- ", std(f)/sqrt(f.n - 1))
            else
                if typeof(f) == Variance{Float64,Float64,EqualWeight}
                    println(k, "\t", mean(f)/avg_sign, " +/- ", std(f)/sqrt(f.n - 1)/avg_sign)
                elseif typeof(f).name == typeof(Group()).name
                    println(string(k))
                    println("-------------")
                    println("means : $(mean.(f.stats) ./ avg_sign)")
                    println("errors : $(std.(f.stats) ./ sqrt(f.n - 1) ./ avg_sign))")
                    println("-------------")
                end
            end
        end
    end

    # calculate kinetic energy in Rydberg
    if in(:Ekin, keys(measurements))
        μT = mean(first(measurements[:Ekin]))/avg_sign
        ΔT = std(first(measurements[:Ekin]))/sqrt(measurements[:Ekin][1].n - 1)/avg_sign
        μT_Ha = E_Ha(μT,e.λ)
        ΔT_Ha = E_Ha(ΔT,e.λ)
        println("T_Ry", "\t", μT_Ha, " +/- ", ΔT_Ha)
    end
    # calculate interaction energy in Rydberg
    if all( [:W_diag,:K_fermion] .∈ (keys(measurements),) )
        μW_diag = mean(first(measurements[:W_diag]))/avg_sign
        ΔW_diag = std(first(measurements[:W_diag]))/sqrt(measurements[:W_diag][1].n - 1)/avg_sign
        μW_off_diag = W_off_diag(e::Ensemble, mean(first(measurements[:K_fermion]))/avg_sign)
        ΔW_off_diag = abs(W_off_diag(e::Ensemble, std(first(measurements[:K_fermion]))/sqrt(measurements[:K_fermion][1].n - 1)/avg_sign))
        μW = μW_diag + μW_off_diag
        ΔW = ΔW_diag + ΔW_off_diag
        ΔWt_Ha = Et_Ha(ΔW, e::Ensemble)
        println("W_diag", "\t", μW_diag, " +/- ", ΔW_diag)
        println("W_off_diag", "\t", μW_off_diag, " +/- ", ΔW_off_diag)
        println("W", "\t", μW, " +/- ", ΔW)
        println("W_t_Ha", "\t", Et_Ha(μW,e), " +/- ", Et_Ha(ΔW,e))
    end
    # calculate total energy in Rydberg
    if all( [:Ekin,:W_diag,:K_fermion] .∈ (keys(measurements),) )
        μE = μW + μT
        ΔE = ΔW + ΔT
        println("E", "\t", μE, " +/- ", ΔE)
    end
end

# TODO: use sign
function save_results(path, measurements, ensemble; options...)
    Nsamples = measurements[:Ekin][1].n

    # thermodynamic observables (single observations)
    fname = joinpath(path, "Obs_N$(ensemble.N)_th$(replace(string(β(ensemble.β,ensemble.N)),"." => ""))_λ$(replace(string(ensemble.ensemble.λ),"." => "")).dat")
    if !isfile(fname)
        print("Print single observables to file ... ")
        open(fname, "a") do io
            writedlm(io, zip(["obs"], ["mean"], ["error"]); options...)
            for (k,(f,m)) in measurements
                if typeof(f) == Variance{Float64,Float64,EqualWeight}
                    # writedlm(io, [typeof(m).name.mt.name, mean(f), std(f) / Nsamples])
                    writedlm(io, zip([string(k)], [mean(f)], [std(f) / Nsamples]); options...)
                end
            end
        end
        print(" . done.\n")
    end

    # structural quantities (group observations)
    for (k,(f,m)) in measurements
        if typeof(f).name == typeof(Group()).name
            fname = joinpath(path, "$(k)_N$(ensemble.N)_th$(replace(string(get_β_internal(ensemble.β,ensemble.N)),"." => ""))_λ$(replace(string(ensemble.λ),"." => "")).dat")
            if !isfile(fname)
                print("Print $(k) to file ... ")
                open(fname, "w") do io
                    writedlm(io, zip(mean.(f.stats), std.(f.stats) ./ Nsamples))
                end
                print(" . done.\n")
            end
        end
    end
end
