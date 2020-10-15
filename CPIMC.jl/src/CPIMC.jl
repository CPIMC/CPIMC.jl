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
    println("Number of Samples : ", Nsamples)
    println("measurements:")
    println("=============")

    for (k,(f,m)) in measurements
        if typeof(f) == Variance{Float64,Float64,EqualWeight}
            println(typeof(m).name.mt.name, "\t", mean(f), " +/- ", std(f) / Nsamples)
        end
    end

    println("")

    println("occupations:")
    println("============")
    println(mean.(measurements[:occs][1].stats))
    println(std.(measurements[:occs][1].stats) ./ Nsamples)
end

function save_results(path, measurements, ensemble)
    Nsamples = measurements[:Ekin][1].n
    # Print to results file
    open(joinpath(path, "occNums_N$(ensemble.N)_th$(replace(string(get_beta_internal(ensemble.beta,ensemble.N)),"." => ""))_rs$(replace(string(ensemble.rs),"." => "")).dat"), "w") do io
           writedlm(io, zip(mean.(measurements[:occs][1].stats), std.(measurements[:occs][1].stats) ./ Nsamples))
       end
end
