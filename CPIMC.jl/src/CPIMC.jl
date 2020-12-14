mutable struct Update
    update! :: Function
    proposed :: Int
    accepted :: Int
end

" change to Configuration "
mutable struct Step{T <: Union{Orbital,Kink,Configuration,Nothing}}
    add :: T
    drop :: T
end

Step{T}() where T = Step{T}(T(),T())
Step() = Step(nothing,nothing)

import Base.empty!# for method extension

# TODO: is this really an in-place operation?
function empty!(Δ)
    Δ = Step()
end

function empty!(Δ::Step{Orbital})
    empty!(Δ.add.occupations)
    empty!(Δ.drop.occupations)
end

function empty!(Δ::Step{Kink})
    empty!(Δ.add.kinks)
    empty!(Δ.drop.kinks)
end

function empty!(Δ::Step{Configuration})
    empty!(Δ.add.kinks)
    empty!(Δ.add.occupations)
    empty!(Δ.drop.kinks)
    empty!(Δ.drop.occupations)
end

import Base.isempty# for method extension

isempty(Δ::Step) = false
isempty(Δ::Step{Nothing}) = true

weight(Δ::Step, e) = weight(Δ.add, Δ.drop, e)

function weight(add::Nothing, drop::Nothing, e) :: Float64
    1.0
end

function weight(add::Orbital, drop::Orbital, e) :: Float64
    return exp( -e.beta * sum( vcat( get_energy(add), -get_energy(drop) ) ) )# explicit summation to minimize cancellation
end

function weight(add::Set{T}, drop::Set{T}, e) :: Float64 where {T}
    return prod(weight.(add, drop))
end

" TODO: implement kink-weight "
function weight(add::SortedDict{Float64, Kink}, drop::SortedDict{Float64, Kink}) :: Float64
    throw(ErrorException("Kinks are not weighted !"))
    return 1.0# TODO: kinks are not weighted !!!
end

" TODO: implement kink-weight "
function weight(add::Kink, drop::Kink, e) :: Float64
    throw(ErrorException("Kinks are not weighted !"))
    return 1.0# TODO: kinks are not weighted !!!
end

function weight(add::Configuration{T}, drop::Configuration{T}, e) :: Float64 where {T}
    return weight(add.occupations, drop.occupations, e) * weight(add.kinks, drop.kinks, e)
end

# " perform multiple MC steps on configuration c "
# function step!(c::Configuration, Δ::Step, e, updates; chain_length)
#     @assert !isempty(updates)
#
#     empty!(Δ)
#
#     δ = Step{basis(c)}()
#
#     chain = zeros(UInt, chain_length)
#     for i in 1:chain_length
#         stp = rand(1:length(updates))
#         chain[i] = stp
#         updates[stp].update!(union(c, Δ), δ, e)
#         Δ.ν *= δ.ν
#         union!(Δ, δ)
#     end
#
#     if rand() < weight(Δ,e) * Δ.ν
#         promote!(c, Δ)
#         for stp in chain
#             updates[stp].proposed += 1
#             updates[stp].accepted += 1
#         end
#     end
# end

" perform an MC step on the configuration c "
function step!(c::Configuration, e::Ensemble, updates::Array{Update,1})
    @assert !isempty(updates)

    up = rand(updates)
    up.proposed += 1
    dv, Δ = up.update!(c, e)
    if rand() < dv
        promote!(c, Δ)
        up.accepted += 1
    end
end

"""change configuration as given by Δ
    TODO: doesn't work with kinks"""
function promote!(c::Configuration, Δ::Step{T}) where T <: Union{Orbital,Kink,Configuration}
    diff!(c, Δ.drop)
    union!(c, Δ.add)
    nothing
end

" empty Step: do nothing "
function promote!(c::Configuration, Δ::Step{Nothing})
    nothing
end

function sweep(steps::Int, sampleEvery::Int, throwAway::Int, updates, measurements, ensemble, c::Configuration; kwargs...)

    " equilibration "
    for i in 1:throwAway
        step!(c, ensemble, updates; kwargs...)
    end

    i = 0

    while i < steps

        " MC step "
        step!(c, ensemble, updates; kwargs...)

        "measurement"
        if i % sampleEvery == 0
            " calculate observables "
            for (key,(stat,obs)) in measurements
                fit!(stat, obs, c)
            end
        end

        i += 1
    end
end

import OnlineStats: fit!# for method extension
fit!(stat::T, obs, c) where T<:Group = fit!(stat, eachrow(obs(c)))
fit!(stat::T, obs, c) where T<:OnlineStat = fit!(stat, obs(c))
