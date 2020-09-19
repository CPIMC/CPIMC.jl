" propose a random update and accept or reject it "
function propose_update!(c::Configuration, updates, e::Ensemble)
    @assert !iszero(length(updates))
    up = rand(updates)
    c_old = Configuration(copy(c.occupations))
    if rand() < up(c,e)
        return :accept
    else
        c.occupations = c_old.occupations
        return :reject
    end
end

" propose a chain of updates "
function propose_update_chain!(c::Configuration, updates, chain_length::Int, e::Ensemble)
end


function runMC(steps::Int, sampleEvery::Int, throwAway::Int, updates, measurements, e::Ensemble, c::Configuration)
    " equilibration "
    for i in 1:throwAway
        propose_update!(c,updates,e)
    end

    i = 0

    while i < steps
        propose_update!(c,updates,e)

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

function save_results(path, measurements, e::Ensemble)
end
