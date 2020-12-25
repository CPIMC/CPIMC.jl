mutable struct Update
    update! :: Function
    proposed :: Int
    accepted :: Int
end

" change to Configuration "
mutable struct Step{S<:Union{Nothing,Basis,Configuration},T<:Union{Nothing,Basis,Configuration}}
    drop :: S
    add :: T
end

Step{T}() where T = Step{T,T}(T(),T())
Step() = Step(nothing,nothing)

import Base.empty!# for method extension

# TODO: is this really an in-place operation?
function empty!(Δ)
    Δ = Step()
end

function empty!(Δ::Step{Configuration})
    empty!(Δ.add.kinks)
    empty!(Δ.add.occupations)
    empty!(Δ.drop.kinks)
    empty!(Δ.drop.occupations)
end

import Base.isempty# for method extension

isempty(Δ::Step) = false
isempty(Δ::Step{Nothing,Nothing}) = true

" return the weight-change of a MC Step. "
weight(Δ::Step, e::Ensemble) :: Float64 = weight(Δ.add, Δ.drop, e)

" return the weight-change of a MC Step given by the changes in the configuration. "
function weight(add::Configuration{T}, drop::Configuration{T}, e::Ensemble) :: Float64 where {T}
    return weight(add.occupations, drop.occupations, e) * weight(add.kinks, drop.kinks, e)
end

" return the weight-change of a MC Step where the configuration is not changed. "
function weight(add::Nothing, drop::Nothing, e::Ensemble) :: Float64
    1.0
end

" return the weight-change of a single diagonal single-particle matrix element "
function weight(add::Orbital, drop::Orbital, e::Ensemble) :: Float64
    exp( -e.beta * ( get_energy(add) - get_energy(drop) ) )
end

" return the weight-change of a set of diagonal single-particle matrix elements "
function weight(add::Set{T}, drop::Set{T}, e::Ensemble) :: Float64 where {T <: Orbital}
    return exp( -e.beta * ( sum( vcat( get_energy.(add),  - get_energy.(drop) ) ) ) )# explicit summation to minimize cancellation
end

" return the weight-change of a SortedDict of off-diagonal single-particle matrix elements (a.k.a. kinks) "
function weight(add::SortedDict{Float64, Kink{T}}, drop::SortedDict{Float64, Kink{T}}, e) :: Float64 where {T}
    return prod( get_energy.(add) ) / prod( get_energy.(drop) )
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
    # !isnothing(Δ.add) && print("proposed step : $(Δ.drop.vec) → $(Δ.add.vec)\nacceptance prob.: $(dv)\tweight: $(weight(Δ, e))")
    # isnothing(Δ.add) && print("no step proposed: $(Δ.drop) -> $(Δ.add)\nacceptance prob.: $(dv)\tweight: $(weight(Δ, e))")
    if rand() < dv
        promote!(c, Δ)
        up.accepted += 1
        # print("   ✓ accepted.\n")
    else
        # print("   x denied.\n")
    end
    # println("current configuration : $([o.vec for o in c.occupations])")

    # println("type $([typeof(o) for o in c.occupations])")
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
