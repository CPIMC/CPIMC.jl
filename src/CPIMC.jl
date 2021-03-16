" a Monte Carlo update to the configuration with two counters "
mutable struct Update
    update :: Function
    proposed :: UInt
    accepted :: UInt
end

" change to Configuration " # TODO Configuration is a mutable type. it may be more efficient to use an immutable dataType for passing the changes in occupations and excitations
struct Step{S<:Union{Nothing,<:Orbital,Configuration},T<:Union{Nothing,<:Orbital,Configuration}}
    drop :: S
    add :: T
end

" outer constructor method for an empty Step "
Step() = Step(nothing,nothing)

# it is possible to define a general outer constructor method Step(d,a) = Step(Configuration(d),Configuration(a))
# this might be bad style and is possibly more difficult to diagnose/debug
# here only specific cases are defined to include the possibility to pass Set{<:Orbital} directly to Step()
# other use-cases are covered by the default constructors
Step(drop::Set{T}, add) where {T <: Orbital} = Step(Configuration(drop), add)
Step(drop, add::Set{T}) where {T <: Orbital} = Step(drop, Configuration(add))
Step(drop::Set{T}, add::Set{T}) where {T <: Orbital} = Step(Configuration(drop), Configuration(add))

" change configuration c as given by Δ "
function apply_step!(c::Configuration, Δ::Step)
    drop!(c, Δ.drop)
    add!(c, Δ.add)
    nothing
end

" change configuration c as given by a list of subsequent Steps "
function apply_step!(c::Configuration, Δ::Array{Step,1})
    for δ in Δ
        apply_step!(c, δ)
    end
end

" return configuration Δ(c) "
apply_step(c::Configuration, Δ::Step) = add(drop(c, Δ.drop), Δ.add)

" change configuration c as given by a list of subsequent Steps "
apply_step(c::Configuration, Δ::Array{Step,1}) = reduce(apply_step, Δ, init=c)

" empty Step: do nothing "
function apply_step!(c::Configuration, Δ::Step{Nothing,Nothing})
    nothing
end

""" perform multiple MC steps on configuration c
    not fully implemented, but to illustrate
    the possibility to achieve this by
    dispatching on the extra_kwarg chain_length """# TODO write function union!(Δ, δ) which combines two Step objects
function update!(c::Configuration, e::Ensemble, updates::Array{Update,1}; chain_length)
    @assert !isempty(updates)

    Δ = Step()
    p = 1

    chain = zeros(UInt, chain_length)
    for i in 1:chain_length
        stp = rand(1:length(updates))
        chain[i] = stp
        dv, δ = updates[stp].update(apply_step(c, Δ), e)
        #TODO write function union!(Δ, δ) which combines two Step objects
        p *= dv
    end

    if rand() < dv
        apply_step!(c, Δ)
        for stp in chain
            updates[stp].proposed += 1
            updates[stp].accepted += 1
        end
    end
end

" perform an MC step on the configuration c "
function update!(c::Configuration, e::Ensemble, updates::Array{Update,1})
    @assert !isempty(updates)

    up = rand(updates)
    up.proposed += 1
    dv, Δ = up.update(c, e)
    if rand() < dv
        apply_step!(c, Δ)
        up.accepted += 1
    end
end


function sweep!(steps::Int, sampleEvery::Int, throwAway::Int, updates::Array{Update,1}, measurements, e::Ensemble, c::Configuration; kwargs...)

    println("starting equilibration")
    k = 1 # progress counter
    for i in 1:throwAway
        # print progress
        if i%(throwAway/100) == 0
            print("eq: ",k,"/100","    ")
            k+=1
        end
        update!(c, e, updates; kwargs...)
    end

    println("starting Simulation")
    i = 0
    k = 1 # progress counter
    while i < steps

        # print progress
        if i%(steps/100) == 0
            print(k,"/100","    ")
            k+=1
        end

        " MC step "
        update!(c, e, updates; kwargs...)

        if i % sampleEvery == 0
            " calculate observables "
            for (key,(stat,obs)) in measurements
                if in(key,[:sign, :K])
                    fit!(stat, obs(e,c))
                else
                    fit!(stat, obs(e,c)*signum(e,c))
                end
            end
        end

        i += 1
    end
end


# TODO: remove
# Only use if all threads do not access any objekts in common
function sweep_multithreaded!(steps::Int, sampleEvery::Int, throwAway::Int, updates::Array{Update,1}, measurements, e::Ensemble, c::Configuration; kwargs...)
    @assert(length(c.kinks) == 0)
    c = Configuration(copy(c.occupations))# c should be a different object for each thread
    " equilibration "
    if (Threads.threadid() == 1)
        println("\nstarting equilibration")
    end
    k = 1# progress counter
    for i in 1:throwAway
        if (i%(throwAway/100) == 0) #& (Threads.threadid() == 1)
            println("               "^(Threads.threadid()-1),"T",Threads.threadid(), " eq: ",k,"/100"," ","K: ",length(c.kinks))
            k+=1
        end
        update!(c, e, updates; kwargs...)
    end
    if (Threads.threadid() == 1)
        println("\nstarting Simulation")
    end
    i = 0
    k = 1#print progress
    #global add_E_counter = 0
    #global remove_E_counter = 0
    while i < steps
        #print progress
        if (i%(steps/100) == 0) #& (Threads.threadid() == 1)
            println("               "^(Threads.threadid()-1),"T",Threads.threadid(), " ",k,"/100"," ","K: ",length(c.kinks))
            k+=1
        end

        " MC step "
        update!(c, e, updates; kwargs...)

        "measurement"
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
                        fit!(stat, obs(e,c)*signum(e,c))
                    end
                end
            end
        end
        i += 1
    end
    println("\nThread",Threads.threadid(),"finished")
end
