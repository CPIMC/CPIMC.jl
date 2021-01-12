" a Monte Carlo update to the configuration with two counters "
mutable struct Update
    update :: Function
    proposed :: UInt
    accepted :: UInt
end

" change to Configuration "# TODO Configuration is a mutable type. it may be more efficient to use an immutable dataType for passing the changes in occupations and excitations
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

""" perform multiple MC steps on configuration c
    not fully implemented, but to illustrate
    the possibility to achieve this by
    dispatching on the extra_kwarg chain_length """# TODO write function union!(Δ, δ) which combines two Step objects
function step!(c::Configuration, e::Ensemble, updates::Array{Update,1}; chain_length)
    @assert !isempty(updates)

    Δ = Step{basis(c)}()
    p = 1

    chain = zeros(UInt, chain_length)
    for i in 1:chain_length
        stp = rand(1:length(updates))
        chain[i] = stp
        dv, δ = updates[stp].update(promote(c, Δ), e)
        #TODO write function union!(Δ, δ) which combines two Step objects
        p *= dv
    end

    if rand() < dv
        promote!(c, Δ)
        for stp in chain
            updates[stp].proposed += 1
            updates[stp].accepted += 1
        end
    end
end

" perform an MC step on the configuration c "
function step!(c::Configuration, e::Ensemble, updates::Array{Update,1})
    @assert !isempty(updates)

    up = rand(updates)
    up.proposed += 1
    dv, Δ = up.update(c, e)
    # println("update :$(up.update)")# debug information
    # !isnothing(Δ.add) && print("proposed step : $(Δ.drop) → $(Δ.add)\nacceptance prob.: $(dv)")# debug information
    # isnothing(Δ.add) && print("no step proposed: $(Δ.drop) → $(Δ.add)\nacceptance prob.: $(dv)")# debug information
    if rand() < dv
        promote!(c, Δ)
        up.accepted += 1
        # print("   ✓ accepted.\n")# debug information
    else
        # print("   x denied.\n")# debug information
    end
    # println("current configuration : $([o.vec for o in c.occupations])")# debug information, TODO: implement printing of configuration in UnicodePlots
end

" change configuration c as given by Δ "
function promote!(c::Configuration, Δ::Step)
    drop!(c, Δ.drop)
    add!(c, Δ.add)
    nothing
end

" return configuration Δ(c) "
promote(c::Configuration, Δ::Step) = add(drop(c, Δ.drop), Δ.add)

" empty Step: do nothing "
function promote!(c::Configuration, Δ::Step{Nothing,Nothing})
    nothing
end

function sweep!(steps::Int, sampleEvery::Int, throwAway::Int, updates::Array{Update,1}, measurements, e::Ensemble, c::Configuration; kwargs...)

    println("starting equilibration")
    k = 1# progress counter
    for i in 1:throwAway
        # print progress
        if i%(throwAway/100) == 0
            print("eq: ",k,"/100","    ")
            k+=1
        end
        step!(c, e, updates; kwargs...)
    end

    println("starting Simulation")
    i = 0
    k = 1# progress counter
    while i < steps

        # print progress
        if i%(steps/100) == 0
            print(k,"/100","    ")
            k+=1
        end

        " MC step "
        step!(c, e, updates; kwargs...)

        if i % sampleEvery == 0
            " calculate observables "
            for (key,(stat,obs)) in measurements
                if typeof(stat) == Group#####################Diese BEdingung ist anscheinend niemals erfüllt
                    fit!(stat, eachrow(obs(e,c)))
                else
                    fit!(stat, obs(e,c))
                end
            end
        end

        i += 1
    end
end


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
        if (i%(throwAway/100) == 0) & (Threads.threadid() == 1)
            print("eq: ",k,"/100","    ")
            k+=1
        end
        step!(c, e, updates; kwargs...)
    end
    if (Threads.threadid() == 1)
        println("\nstarting Simulation")
    end
    i = 0
    k = 1# progress counter
    while i < steps
        # print progress
        if (i%(steps/100) == 0) & (Threads.threadid() == 1)
            print(k,"/100","    ")
            #println("K: ",length(c.kinks))# debug information
            k+=1
        end

        " MC step "
        step!(c, e, updates; kwargs...)

        "measurement"
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
    println("\nThread",Threads.threadid(),"finished")
end
