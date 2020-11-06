" propose a random update and accept or reject it "
function propose_update!(c::Configuration, updates, ensemble)
    @assert !iszero(length(updates))
    up = rand(updates)
    c_old = Configuration(copy(c.occupations))
    if rand() < up(c,ensemble)
        return :accept
    else
        c.occupations = c_old.occupations
        return :reject
    end
end

" propose a chain of updates "
function propose_update_chain!(c::Configuration, updates, chain_length::Int, ensemble)
end


function runMC(steps::Int, sampleEvery::Int, throwAway::Int, updates, measurements, ensemble, c::Configuration)
    " equilibration "
    for i in 1:throwAway
        propose_update!(c,updates,ensemble)
    end

    i = 0

    while i < steps
        propose_update!(c,updates,ensemble)

        if i % sampleEvery == 0
            " calculate observables "
            for (key,(stat,obs)) in measurements
                if typeof(stat) == Group
                    fit!(stat, eachrow(obs(c)))
                else
                    fit!(stat, obs(c))
                end
            end
        end

        i += 1
    end
end

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
