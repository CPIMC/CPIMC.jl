" propose a random update and accept or reject it "
function propose_update!(c::Configuration, updates, e::Ensemble)
    @assert !iszero(length(updates))
    up = rand(updates)
    #print(up)
    c_old = Configuration(copy(c.occupations),copy(c.kinks))
    if rand() < up(c,e)
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
function propose_update_chain!(c::Configuration, updates, chain_length::Int, e::Ensemble)
end


function runMC(steps::Int, sampleEvery::Int, throwAway::Int, updates, measurements, e::Ensemble, c::Configuration)
    " equilibration "
    println("starting equilibration")
    k = 1 #print prgoress
    for i in 1:throwAway
        if i%(throwAway/100) == 0
            println("eq: ",k,"/100")
            k+=1
        end
        propose_update!(c,updates,e)
    end
    println("starting Simulation")
    i = 0
    k = 1#print progress
    while i < steps
        #print progress
        if i%(steps/100) == 0
            println(k,"/100")
            k+=1
        end
        propose_update!(c,updates,e)
        """try###Debugg Code
            get_occupations_at(c, e.beta)
        catch a
            "brakepoint"
        end"""
        if i % sampleEvery == 0
            " calculate observables "
            for (key,(stat,obs)) in measurements
                if typeof(stat) == Group
                    fit!(stat, eachrow(obs(e,c)))
                else
                    fit!(stat, obs(e,c))
                end
            end
        end

        i += 1
    end
end

function save_results(path, measurements, e)
end
