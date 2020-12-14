
function print_results(measurements)
    Nsamples = measurements[:Ekin][1].n
    println("number of samples : ", Nsamples)
    println("measurements:")
    println("=============")

    for (k,(f,m)) in measurements
        if typeof(f) == Variance{Float64,Float64,EqualWeight}
            println(string(k), "\t", mean(f), " +/- ", std(f) / Nsamples)
        elseif typeof(f).name == typeof(Group()).name
            println(string(k))
            println("-------------")
            println("means : $(mean.(f.stats))")
            println("errors : $(std.(f.stats) ./ Nsamples))")
            println("-------------")
        end
    end
end

function save_results(path, measurements, ensemble; options...)
    Nsamples = measurements[:Ekin][1].n

    # thermodynamic observables (single observations)
    fname = joinpath(path, "Obs_N$(ensemble.N)_th$(replace(string(get_beta_internal(ensemble.beta,ensemble.N)),"." => ""))_rs$(replace(string(ensemble.rs),"." => "")).dat")
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
            fname = joinpath(path, "$(k)_N$(ensemble.N)_th$(replace(string(get_beta_internal(ensemble.beta,ensemble.N)),"." => ""))_rs$(replace(string(ensemble.rs),"." => "")).dat")
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
