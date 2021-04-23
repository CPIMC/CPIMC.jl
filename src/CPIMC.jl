
"""
Common functionality to simulate any physical model.
"""
module CPIMC
using StaticArrays
using DataStructures
using FixedPointNumbers
import LinearAlgebra: dot, norm
using OnlineStats
using Printf

export Group, Mean, Variance, Fixed


include("configuration.jl")
include("ensemble.jl")
include("model.jl")

include("Estimators.jl")
include("PlaneWaves.jl")
include("DefaultUpdates.jl")
include("UniformElectronGas.jl")


export UpdateCounter, sweep!, print_results, print_rates, Step, apply_step, apply_step!, apply_back_step!, update!, measure!

"""
Objects of this type are used to keep track of the acceptance ratio of a class of updates.

**Fields**
- `rejected` -- number of times updates have been rejected
- `accepted` -- number of times updates have been accepted
- `trivial`  -- number of times the changes proposed by updates of that class have been trivial, e.g. the Step returned by update has been empty

For convenience, an outer constructor with all counters set to zero is provided:

    UpdateCounter() = UpdateCounter(0,0,0)

"""
mutable struct UpdateCounter
    rejected :: Int64
    accepted :: Int64
    trivial  :: Int64
end

UpdateCounter() = UpdateCounter(0,0,0)

"""
    Base.:+(x::UpdateCounter, y::UpdateCounter)

Addition of counters for example to add counters from different Markov-chains.
"""
function Base.:+(x::UpdateCounter, y::UpdateCounter)
    UpdateCounter(x.rejected + y.rejected, x.accepted + y.accepted, x.trivial + y.trivial)
end

"""
    Diff = Union{Nothing,Tuple}

Type alias for Union{Nothing,Tuple}. This is used as type parameters for the fields of struct Step.
Mainly for convenient construction of any combination of ::Nothing, ::Tuple.
"""
const Diff = Union{Nothing,Tuple}
Diff(::Nothing) = nothing

"""
Represents a change which can be applied to a `Configuration`.

**Fields**

- `drop_orbs`  -- `Orbital`s which are to be removed from the initial occupation
- `drop_kinks` -- pairs of `ImgTime` and `Kink`s which are to be removed
- `add_orbs`   -- `Orbital`s which are added to the initial occupation
- `add_kinks`  -- pairs of `ImgTime` and `Kinks` which are added

"""
struct Step{T<:Diff, S<:Diff, R<:Diff, Q<:Diff}
    drop_orbs :: T
    drop_kinks :: S
    add_orbs :: R
    add_kinks :: Q
end


# outer constructor method for an empty Step
Step() = Step(nothing,nothing,nothing,nothing)
Step(a,b) = Step(a,nothing,b,nothing)

"""
    apply_step!(c::Configuration, step::Step)
    apply_step!(c::Configuration, steps)

Apply the changes given by the second argument to a `Configuration`.
"""
function apply_step!(c::Configuration, s::Step)
    drop_orbs!(c.occupations, s.drop_orbs)
    drop_kinks!(c.kinks, s.drop_kinks)
    add_orbs!(c.occupations, s.add_orbs)
    add_kinks!(c.kinks, s.add_kinks)
end

function apply_back_step!(c::Configuration, s::Step)
    drop_orbs!(c.occupations, s.add_orbs)
    drop_kinks!(c.kinks, s.add_kinks)
    add_orbs!(c.occupations, s.drop_orbs)
    add_kinks!(c.kinks, s.drop_kinks)
end

function apply_step!(c::Configuration, steps::Array)
    for step in steps
        apply_step!(c, step)
    end
end

function apply_back_step!(c::Configuration, steps::Array)
    for i in 0:-1:-(length(steps)-1)
        apply_back_step!(c, steps[length(steps)-i])
    end
end

"""
    apply_step(c::Configuration, step::Step)
    apply_step(c::Configuration, steps)

Return the result of applying the changes given by the second argument to a `Configuration`.
"""
function apply_step(c::Configuration, steps)
    c1 = deepcopy(c)
    apply_step!(c1,steps)
    return c1
end


"""
    update!(m::Model, e::Ensemble, c::Configuration, updates)

Obtain the next state of the Markov chain by randomly selecting an element of `updates`
and applying the proposed changes to `c` with the calculated probability. Returns a tuple
containing the index of the chosen function and `:accept`, `:reject` or `:trivial` if
the function yielded an empty `Step`.
"""
function update!(m::Model, e::Ensemble, c::Configuration, updates)
    @assert !isempty(updates)
    i = rand(1:length(updates))
    dv, Δ = updates[i](m, e, c)

    if Δ == Step()
        return (i, :trivial)
    end

    if rand() < dv
        return (i, :accept)
    else
        apply_back_step!(c, Δ)
        return (i, :reject)
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
    sweep!(m::Model, e::Ensemble, c::Configuration, updates, estimators, steps::Int, sampleEvery::Int, throwAway::Int)

Generate a markov chain of length `steps` using the Metropolis-Hastings algorithm with the updates given in `updates`.
After `throwAway` steps have been performed, the observables given in `estimators` are calculated every `sampleEvery` steps.

Returns a dictionary containing an `UpdateCounter` for every function given in `updates`.
"""
function sweep!(m::Model, e::Ensemble, c::Configuration, updates, estimators, steps::Int, sampleEvery::Int, throwAway::Int)
    counters = [ UpdateCounter() for u in updates ]

    if (Threads.threadid() == 1)
        println("\nstarting equilibration")
    end

    k = 1 # progress counter
    for i in 1:throwAway
        # print progress
        if i % (throwAway/100) == 0
            print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀    "^(Threads.threadid()-1), "eq:",Threads.threadid(), "T",lpad(k,3, ' '),"%","K:",lpad(length(c.kinks),4, ' '), "\r")
            k += 1
        end
        update!(m, e, c, updates)
    end
    if (Threads.threadid() == 1)
        println("\nstarting simulation")
    end

    i = 0
    k = 1 # progress counter

    while i < steps
        # print progress
        if i % (steps/100) == 0
            print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀    "^(Threads.threadid()-1),"   ",Threads.threadid(), "T",lpad(k,3, ' '),"%","K:",lpad(length(c.kinks),4, ' '), "\r")
            k+=1
        end

        " MC step "
        (up, res) = update!(m, e, c, updates)

        if res == :accept
            counters[up].accepted += 1
        elseif res == :reject
            counters[up].rejected += 1
        elseif res == :trivial
            counters[up].trivial += 1
        end

        " calculate estimators "
        if i % sampleEvery == 0
            measure!(m, e, c, estimators)
        end

        i += 1
    end
    println("\nthread ",Threads.threadid()," finished")

    return Dict(zip(map(u -> typeof(u).name.mt.name, updates), counters))
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

"""
    print_rates(dict::Dict{Any,UpdateCounter})

Calculate and print acceptance statistics.
"""
function print_rates(dict)
    long = 31

    println("\n")
    println(rpad("update", long), lpad("accepted", 11), lpad("rejected", 11), lpad("trivial", 11))
    println("================================================================")
    println("\ntotal:\n")

    for (up,cn) in dict

        println(rpad(up, long), lpad(cn.accepted, 11), lpad(cn.rejected, 11), lpad(cn.trivial, 11))
    end
    overall_counters = sum(values(dict))
    println(rpad("all", long), lpad(overall_counters.accepted, 11), lpad(overall_counters.rejected, 11), lpad(overall_counters.trivial, 11))
    println("\nrelative:\n")

    for (up,cn) in dict
        proposed = cn.accepted + cn.rejected + cn.trivial
        println(rpad(up, long),
                lpad(@sprintf("%.2f", 100*cn.accepted/proposed), 10), "%",
                lpad(@sprintf("%.2f", 100*cn.rejected/proposed), 10), "%",
                lpad(@sprintf("%.2f", 100*cn.trivial/proposed), 10), "%")
    end
    proposed = overall_counters.accepted + overall_counters.rejected + overall_counters.trivial
    println(rpad("all", long),
            lpad(@sprintf("%.2f", 100*overall_counters.accepted/proposed), 10), "%",
            lpad(@sprintf("%.2f", 100*overall_counters.rejected/proposed), 10), "%",
            lpad(@sprintf("%.2f", 100*overall_counters.trivial/proposed), 10), "%")
    println("\n\n")
end

end
