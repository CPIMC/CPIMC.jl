"""
new implementation of the CPIMC method
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


" a Monte Carlo update to the configuration with two counters "
mutable struct Update
    update :: Function
    proposed :: UInt
    accepted :: UInt
    trivial :: UInt
end

Update(f::Function) = Update(f,0,0,0)


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
apply_step(c::Configuration, Δ::Array{Step,1}) = foldl(apply_step, Δ, init=c)

" empty Step: do nothing "
function apply_step!(c::Configuration, Δ::Step{Nothing,Nothing})
    nothing
end

" perform a MC step on the configuration c "
function update!(m::Model, e::Ensemble, c::Configuration, updates::Array{Update,1})
    @assert !isempty(updates)
    up = rand(updates)
    up.proposed += 1
    dv, Δ = up.update(m, e, c)
    if rand() < dv
        apply_step!(c, Δ)
        if Δ == Step()
            up.trivial += 1
        else
            up.accepted += 1
        end
    end
end

"""
    measure!(m::Model, e::Ensemble, c::Configuration, estimators)

Perform measurements on a configuration. `estimators` needs to be a dictionary that contains tuples `(::OnlineStat, ::Function)`.
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

    println("starting equilibration")
    k = 1 # progress counter
    for i in 1:throwAway
        # print progress
        if i % (throwAway/100) == 0
            print("eq: ", k, "/100   ", "K ", lpad(length(c.kinks),4, ' '), "\r")
            k += 1
        end
        update!(m, e, c, updates)
    end

    println("\n starting Simulation")
    i = 0
    k = 1 # progress counter
    while i < steps
        # print progress
        if i % (steps/100) == 0
            print(k, "/100   ", "K: ", lpad(length(c.kinks),4, ' '), "\r")
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

    println("\n")
end

" print all measured observables and calculate further quantities "
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
