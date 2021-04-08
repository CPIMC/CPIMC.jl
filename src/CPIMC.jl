"""
Common functionality to simulate any physical model.
"""
module CPIMC

using DataStructures
using FixedPointNumbers
import LinearAlgebra: dot, norm
using OnlineStats

export Group, Mean, Variance


include("configuration.jl")
include("ensemble.jl")
include("model.jl")

include("Estimators.jl")
include("PlaneWaves.jl")


"""
    const ex_radius = 3

radius of the sphere of orbitals which are considered for excitations
"""
const ex_radius = 3 # TODO: find better solution
# ex_radius could be kwarg to each update function
# and passed to `sweep!` via anonymous function with local variable ex_radius


using .PlaneWaves

include("updates/type_a.jl")
include("updates/type_b.jl")
include("updates/type_c.jl")
include("updates/type_d.jl")
include("updates/type_e.jl")

include("updates/auxiliary.jl")


include("UniformElectronGas.jl")


export Update, sweep!, print_results


"""
Wrapper structure for a Monte Carlo update.

**Fields**

- `update :: Function` -- actual update function
- `proposed`           -- number of times this update has been proposed
- `accepted`           -- number of times this update has been accepted
- `trivial`            -- number of times the changes proposed by this update have been trivial, e.g. the `Step` returned by `update` has been empty

For convenience, an outer constructor with all counters set to zero is provided:

    Update(f::Function) = Update(f,0,0,0)

"""
mutable struct Update  
    update :: Function
    proposed :: UInt
    accepted :: UInt
    trivial :: UInt
end

Update(f::Function) = Update(f,0,0,0)


"""
Represents a change which can be applied to a `Configuration`.

**Fields**

- `drop` 
- `add`

"""
struct Step{S<:Union{Nothing,<:Orbital,Configuration},T<:Union{Nothing,<:Orbital,Configuration}}
    drop :: S
    add :: T
end
# TODO Configuration is a mutable type. it may be more efficient to use an immutable dataType for passing the changes in occupations and excitations

# outer constructor method for an empty Step
Step() = Step(nothing,nothing)


# it is possible to define a general outer constructor method Step(d,a) = Step(Configuration(d),Configuration(a))
# this might be bad style and is possibly more difficult to diagnose/debug
# here only specific cases are defined to include the possibility to pass Set{<:Orbital} directly to Step()
# other use-cases are covered by the default constructors
Step(drop::Set{T}, add) where {T <: Orbital} = Step(Configuration(drop), add)
Step(drop, add::Set{T}) where {T <: Orbital} = Step(drop, Configuration(add))
Step(drop::Set{T}, add::Set{T}) where {T <: Orbital} = Step(Configuration(drop), Configuration(add))

"""
    apply_step!(c::Configuration, step::Step)
    apply_step!(c::Configuration, steps)

Apply the changes given by the second argument to a `Configuration`.
"""
function apply_step!(c::Configuration, step::Step)
    drop!(c, step.drop)
    add!(c, step.add)
    nothing
end

function apply_step!(c::Configuration, steps)
    for step in steps
        apply_step!(c, step)
    end
end

"""
    apply_step(c::Configuration, step::Step)
    apply_step(c::Configuration, steps)

Return the result of applying the changes given by the second argument to a `Configuration`.
"""
apply_step(c::Configuration, step::Step) = add(drop(c, step.drop), step.add)
apply_step(c::Configuration, steps) = foldl(apply_step, steps, init=c)

#empty Step: do nothing
function apply_step!(c::Configuration, Δ::Step{Nothing,Nothing})
    nothing
end

"""
    update!(m::Model, e::Ensemble, c::Configuration, updates::Array{Update,1})

Obtain the next state of the Markov chain by randomly selecting an element of `updates`
and applying the proposed changes to `c` with the calculated probability.
"""
function update!(m::Model, e::Ensemble, c::Configuration, updates::Array{Update,1})
    @assert !isempty(updates)
    no_kinks_Updates = filter(x -> in(x.update,[move_particle,add_type_B]), updates)
    two_kinks_Updates = filter(x -> in(x.update,[move_particle, add_type_B, remove_type_B, add_type_C, add_type_D, add_type_E, add_remove_kink_chain, shuffle_indices]), updates)
    threeplus_kinks_Updates = filter(x -> in(x.update,[move_particle, add_type_B, remove_type_B, add_type_C, remove_type_C, add_type_D, remove_type_D, add_type_E, remove_type_E, add_remove_kink_chain, shuffle_indices]), updates)
    if isempty(c.kinks)
        up = rand(no_kinks_Updates)
        up.proposed += 1
        dv, Δ = up.update(m, e, c)
        if !isempty(apply_step(c, Δ).kinks)
            @assert(length(apply_step(c, Δ).kinks) == 2)
            dv *= length(no_kinks_Updates)/length(two_kinks_Updates)
        end
        if rand() < dv
            apply_step!(c, Δ)
            if Δ == Step()
                up.trivial += 1
            else
                up.accepted += 1
            end
        end
    elseif length(c.kinks) == 2
        up = rand(two_kinks_Updates)
        up.proposed += 1
        dv, Δ = up.update(m, e, c)
        if isempty(apply_step(c, Δ).kinks)
            dv *= length(two_kinks_Updates)/length(no_kinks_Updates)
        elseif length(apply_step(c, Δ).kinks) != 2
            dv *= length(two_kinks_Updates)/length(threeplus_kinks_Updates)
        end
        if rand() < dv
            apply_step!(c, Δ)
            if Δ == Step()
                up.trivial += 1
            else
                up.accepted += 1
            end
        end
    else
        up = rand(threeplus_kinks_Updates)
        up.proposed += 1
        dv, Δ = up.update(m, e, c)
        if length(apply_step(c, Δ).kinks) == 2
            dv *= length(threeplus_kinks_Updates)/length(two_kinks_Updates)
        end
        if rand() < dv
            apply_step!(c, Δ)
            if Δ == Step()
                up.trivial += 1
            else
                up.accepted += 1
            end
        end
    end
end

"""
    measure!(m::Model, e::Ensemble, c::Configuration, estimators)

Perform measurements on a `Configuration`. `estimators` needs to be a dictionary that contains tuples `(::OnlineStat, ::Function)`.
"""
function measure!(m, e, c, estimators)
    s = signum(m, c)
    for (key, (stat, obs)) in estimators
        if in(key, [:sign, :K, :T1c])
            fit!(stat, obs(m,e,c))
        else
            fit!(stat, obs(m,e,c)*s)
        end
    end
end


"""
    sweep!(m::Model, e::Ensemble, c::Configuration, updates::Array{Update,1}, measurements, steps::Int, sampleEvery::Int, throwAway::Int)

Generate a markov chain of length `steps` using the Metropolis-Hastings algorithm with the updates given in `updates`.
After `throwAway` steps have been performed, the observables given in `estimators` are calculated every `sampleEvery` steps.
"""
function sweep!(m::Model, e::Ensemble, c::Configuration, updates::Array{Update,1}, estimators, steps::Int, sampleEvery::Int, throwAway::Int)

    if (Threads.threadid() == 1)
        println("\nstarting equilibration")
    end
    k = 1 # progress counter
    for i in 1:throwAway
        # print progress
        if i % (throwAway/100) == 0
            print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀    "^(Threads.threadid()-1), "eq:",Threads.threadid(), "T",k,"%"," K:",lpad(length(c.kinks),4, ' '), "\r")
            k += 1
        end
        update!(m, e, c, updates)
    end
    if (Threads.threadid() == 1)
        println("\n starting Simulation")
    end
    i = 0
    k = 1 # progress counter
    while i < steps
        # print progress
        if i % (steps/100) == 0
            print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀ ⠀   "^(Threads.threadid()-1),"   ",Threads.threadid(), "T",k,"%"," K:",lpad(length(c.kinks),4, ' '), "\r")
            k+=1
        end

        " MC step "
        update!(m, e, c, updates)

        " calculate estimators "
        if i % sampleEvery == 0
            measure!(m, e, c, estimators)
        end

        i += 1
    end
    println("\nThread",Threads.threadid(),"finished")
end

"""
    print_results(measurements, e::Ensemble)

Print all measured observables, potentially taking the sign into account.
"""
function print_results(measurements, e::Ensemble)
    println("measurements:")
    println("=============")

    # calculate average sign
    if haskey(measurements, :sign)
        avg_sign = mean(first(measurements[:sign]))
    else
        avg_sign = 1.0
    end

    for (k,(f,m)) in measurements
        if typeof(f) == Variance{Float64,Float64,EqualWeight}
            if in(k,[:sign, :K])
                println(k, "\t", mean(f), " +/- ", std(f)/sqrt(f.n - 1))
            else
                if typeof(f) == Variance{Float64,Float64,EqualWeight}
                    println(k, "\t", mean(f)/avg_sign, " +/- ", std(f)/sqrt(f.n - 1)/avg_sign)
                elseif typeof(f).name == typeof(Group()).name
                    println(string(k))
                    println("-------------")
                    println("means : $(mean.(f.stats) ./ avg_sign)")
                    println("errors : $(std.(f.stats) ./ sqrt(f.n - 1) ./ avg_sign))")
                    println("-------------")
                end
            end
        end
    end
end


end
