" propose a random update and accept or reject it "
function propose_update!(c::Configuration, updates, e::Ensemble)
    @assert !iszero(length(updates))
    up = rand(updates)
    c_old = Configuration(copy(c.occupations),copy(c.kinks))
    if rand() < up(c,e)
        return :accept
    else
        c.occupations = c_old.occupations
        c.kinks = c_old.kinks
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
            print("eq: ",k,"/100","    ")
            #println("K: ",length(c.kinks))
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
            print(k,"/100","    ")
            #println("K: ",length(c.kinks))
            k+=1
        end
        propose_update!(c,updates,e)
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

#measurements should be the only variable modified by multyple threads
function runMC_multythreaded(steps::Int, sampleEvery::Int, throwAway::Int, updates, measurements, e::Ensemble, c::Configuration)
    l = Threads.ReentrantLock()
    @assert(length(c.kinks) == 0)
    c = Configuration(copy(c.occupations))  #c should be a different objekt for each thread
    " equilibration "
    if (Threads.threadid() == 1)
        println("\nstarting equilibration")
    end
    k = 1 #print prgoress
    for i in 1:throwAway
        if (i%(throwAway/100) == 0) & (Threads.threadid() == 1)
            print("eq: ",k,"/100","    ")
            #println("K: ",length(c.kinks))
            k+=1
        end
        propose_update!(c,updates,e)
    end
    if (Threads.threadid() == 1)
        println("\nstarting Simulation")
    end
    i = 0
    k = 1#print progress
    while i < steps
        #print progress
        if (i%(steps/100) == 0) & (Threads.threadid() == 1)
            print(k,"/100","    ")
            #println("K: ",length(c.kinks))
            k+=1
        end
        propose_update!(c,updates,e)
        if i % sampleEvery == 0
            " calculate observables "
            for (key,(stat,obs)) in measurements
                if typeof(stat) == Group
                    Threads.lock(l)
                        fit!(stat, eachrow(obs(e,c)))
                    Threads.unlock(l)
                else
                    Threads.lock(l)
                        fit!(stat, obs(e,c))
                    Threads.unlock(l)
                end
            end
        end
        i += 1
    end
    println("\nThread",Threads.threadid(),"finished")
end


function save_results(path, measurements, e)
end
