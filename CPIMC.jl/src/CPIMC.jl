#global add_c_counter = 0
#global remove_c_counter = 0

" propose a random update and accept or reject it "
function propose_update!(c::Configuration, updates, e::Ensemble)
    @assert !iszero(length(updates))
    up = rand(updates)

    c_old = Configuration(copy(c.occupations),copy(c.kinks),copy(c.sign))
    acc_prob = up(c,e)  #Dies in Variable zu sspeichern nur sinnvoll fürs debuggen
    if (rand() < acc_prob) #& (acc_prob != 1)
        """if up == add_type_C
            global add_c_counter += 1
        end
        if up == remove_type_C
            global remove_c_counter += 1
        end"""
        return :accept
    else
        """if acc_prob == 1
            return :accept
        end"""
        c.occupations = c_old.occupations
        c.kinks = c_old.kinks
        c.sign = c_old.sign
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
            println("K: ",length(c.kinks))
            k+=1
        end
        propose_update!(c,updates,e)
    end
    println("starting Simulation")
    global add_c_counter = 0
    global remove_c_counter = 0
    i = 0
    k = 1#print progress
    while i < steps
        #print progress
        if i%(steps/100) == 0
            print(k,"/100","    ")
            println("K: ",length(c.kinks))
            k+=1
        end
        propose_update!(c,updates,e)
        if i % sampleEvery == 0
            " calculate observables "
            for (key,(stat,obs)) in measurements
                if in(key,[:sign, :K])
                    fit!(stat, obs(e,c))
                else
                    if typeof(stat) == Group#####################Diese Bedingung ist anscheinend niemals erfüllt
                        println("Das wird nicht geprinted")
                        fit!(stat, eachrow(obs(e,c)))
                    else
                        fit!(stat, obs(e,c)*c.sign)
                    end
                end
            end
        end

        i += 1
    end
    #println(add_c_counter)
    #println(remove_c_counter)
end


#Only use if all threads do not access any objekts in common
function runMC_multithreaded(steps::Int, sampleEvery::Int, throwAway::Int, updates, measurements, e::Ensemble, c::Configuration)
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
            println("K: ",length(c.kinks))
            k+=1
        end
        propose_update!(c,updates,e)
    end
    if (Threads.threadid() == 1)
        println("\nstarting Simulation")
    end
    i = 0
    k = 1#print progress
    global add_c_counter = 0
    global remove_c_counter = 0
    while i < steps
        #print progress
        if (i%(steps/100) == 0) & (Threads.threadid() == 1)
            print(k,"/100","    ")
            println("K: ",length(c.kinks))
            k+=1
        end
        propose_update!(c,updates,e)
        if i % sampleEvery == 0
            " calculate observables "
            for (key,(stat,obs)) in measurements
                if in(key,[:sign, :K])
                    fit!(stat, obs(e,c))
                else
                    if typeof(stat) == Group#####################Diese Bedingung ist anscheinend niemals erfüllt
                        println("Das wird nicht geprinted")
                        fit!(stat, eachrow(obs(e,c)))
                    else
                        fit!(stat, obs(e,c)*c.sign)
                    end
                end
            end
        end
        i += 1
    end
    println("\nThread",Threads.threadid(),"finished")
    #println(add_c_counter)
    #println(remove_c_counter)
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
