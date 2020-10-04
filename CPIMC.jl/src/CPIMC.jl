" propose a random update and accept or reject it "
function propose_update!(c::Configuration, updates, ensemble)
    @assert !iszero(length(updates))
    up = rand(updates)
    #print(up)
    c_old = Configuration(copy(c.occupations),copy(c.kinks))
    if rand() < up(c,ensemble)
        #print("   accepted\n")
        return :accept
    else
        c.occupations = c_old.occupations
        c.kinks = c_old.kinks
        #print("   rejected\n")
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
        """try###Debugg Code
            get_occupations_at(c, ensemble.beta)
        catch a
            "brakepoint"
        end"""
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

function save_results(path, measurements, ensemble)
end
