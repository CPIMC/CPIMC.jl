using StaticArrays

"""
    Diff = Union{Nothing,NTuple}

Type alias for Union{Nothing,NTuple}. This is used as type parameters for the fields of struct Step.
Mainly for convenient construction of any combination of ::Nothing, ::TTuple.
"""
const Diff_svec = Union{Nothing,SVector}
Diff_svec(::Nothing) = nothing
Diff_svec(t::Tuple) = SVector(t)

const Diff_tuple = Union{Nothing,Tuple}
Diff_tuple(::Nothing) = nothing

struct StepA{T<:Diff_svec, S<:Diff_svec, R<:Diff_svec, Q<:Diff_svec}
    drop_orbs :: T
    drop_kinks :: S
    add_orbs :: R
    add_kinks :: Q
end

StepA() = StepA(nothing, nothing, nothing, nothing)
"""
    StepA(a,b,c,d)

Construct a Step for any combination of ::Nothing and ::SVector.
"""
# Without using Diff_svec it would be needed to define 4! constructors for all combinations
StepA(a,b,c,d) = StepA(Diff_svec(a), Diff_svec(b), Diff_svec(c), Diff_svec(d))

struct StepB{T<:Diff_tuple, S<:Diff_tuple, R<:Diff_tuple, Q<:Diff_tuple}
    drop_orbs :: T
    drop_kinks :: S
    add_orbs :: R
    add_kinks :: Q
end

StepB() = StepB(nothing, nothing, nothing, nothing)
StepB(a,b,c,d) = StepB(Diff_tuple(a), Diff_tuple(b), Diff_tuple(c), Diff_tuple(d))

struct StepC
    drop_orbs :: Union{SVector, Nothing}
    drop_kinks :: Union{SVector, Nothing}
    add_orbs :: Union{SVector, Nothing}
    add_kinks :: Union{SVector, Nothing}
end

StepC() = StepC(nothing, nothing, nothing, nothing)

struct StepD{T<:SVector, S<:SVector, R<:SVector, Q<:SVector}
    drop_orbs :: T
    drop_kinks :: S
    add_orbs :: R
    add_kinks :: Q
end

emptyStepD(T::DataType) = StepD(SVector{0,T}(),SVector{0,T}(),SVector{0,T}(),SVector{0,T}())
StepD(a,b,c,d) = StepD(Diff_svec(a), Diff_svec(b), Diff_svec(c), Diff_svec(d))

struct StepE
    drop_orbs :: SVector
    drop_kinks :: SVector
    add_orbs :: SVector
    add_kinks :: SVector
end

emptyStepE(T::DataType) = StepE(SVector{0,T}(),SVector{0,T}(),SVector{0,T}(),SVector{0,T}())

struct StepF{T<:Tuple, S<:Tuple, R<:Tuple, Q<:Tuple}
    drop_orbs :: T
    drop_kinks :: S
    add_orbs :: R
    add_kinks :: Q
end

StepF() = StepF((),(),(),())
StepF(a,b,c,d) = StepF(Diff_tuple(a), Diff_tuple(b), Diff_tuple(c), Diff_tuple(d))

struct StepF{T<:Tuple, S<:Tuple, R<:Tuple, Q<:Tuple}
    drop_orbs :: T
    drop_kinks :: S
    add_orbs :: R
    add_kinks :: Q
end

StepF() = StepF((),(),(),())
StepF(a,b,c,d) = StepF(Diff_tuple(a), Diff_tuple(b), Diff_tuple(c), Diff_tuple(d))

struct StepG
    drop_orbs :: Tuple
    drop_kinks :: Tuple
    add_orbs :: Tuple
    add_kinks :: Tuple
end

StepG() = StepG((),(),(),())

drop_orbs!(occ, ::Nothing) = nothing
drop_orbs!(occ, s::SVector{0,T}) where T = nothing
drop_orbs!(occ, s::Tuple{}) where T = nothing
drop_orbs!(occ, s) = setdiff!(occ, s)


add_orbs!(occ, ::Nothing) = nothing
add_orbs!(occ, s::SVector{0,T}) where T = nothing
add_orbs!(occ, s::Tuple{}) where T = nothing
add_orbs!(occ, s) = union!(occ, s)

drop_kinks!(ck, ::Nothing) = nothing
drop_kinks!(ck, s::SVector{0,T}) where T = nothing
drop_kinks!(ck, s::Tuple{}) where T = nothing
drop_kinks!(ck, s) = setdiff!(ck , s)

add_kinks!(ck, ::Nothing) = nothing
add_orbs!(ck, s::SVector{0,T}) where T = nothing
add_orbs!(ck, s::Tuple{}) where T = nothing
function add_kinks!(ck, s)
    for k in s
        insert!(ck, searchsortedfirst(ck, by=first, first(k)), first(k) => last(k))
    end
end

function apply_step!(c, s)
    drop_orbs!(c.occupations, s.drop_orbs)
    drop_kinks!(c.kinks, s.drop_kinks)
    add_orbs!(c.occupations, s.add_orbs)
    add_kinks!(c.kinks, s.add_kinks)
end

import CPIMC: T2, T4, Kink, ImgTime, Configuration, Kinks, Orbital

struct Step{S<:Union{Nothing,<:Orbital,Configuration},T<:Union{Nothing,<:Orbital,Configuration}}
    drop :: S
    add :: T
end

Step() = Step(nothing,nothing)
Step(drop::Set{T}, add) where {T <: Orbital} = Step(Configuration(drop), add)
Step(drop, add::Set{T}) where {T <: Orbital} = Step(drop, Configuration(add))
Step(drop::Set{T}, add::Set{T}) where {T <: Orbital} = Step(Configuration(drop), Configuration(add))

function apply_step!(c::Configuration, step::Step)
    drop!(c, step.drop)
    add!(c, step.add)
    nothing
end


drop!(c::Configuration, n::Nothing) = nothing

drop!(c::Configuration{T}, o::T) where {T <: Orbital} = delete!(c.occupations, o)
function drop!(c::Configuration, oc)
    for o in oc
        drop!(c, o)
    end
end
drop!(c::Configuration, τ::ImgTime) = setdiff!(c.kinks, [τ => c.kinks[τ]])
drop!(c::Configuration, τ, k) = setdiff!(c.kinks, τ => k)
function drop!(c::Configuration{T}, ck::Kinks) where {T <: Orbital}
    for (τ, k) in ck
        drop!(c, τ, k)
    end
end
function drop!(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital}
    drop!(c1, c2.occupations)
    drop!(c1, c2.kinks)
end
function add!(c::Configuration, n::Nothing)
    nothing
end

add!(c::Configuration{T}, o::T) where {T <: Orbital} = push!(c.occupations, o)
add!(c::Configuration, oc) = union!(c.occupations, oc)
add!(ck::Kinks, τ::ImgTime, k::Kink) = insert!(ck, searchsortedfirst(ck, by=first, τ), τ => k)
add!(c::Configuration{T}, p::Pair ) where {T <: Orbital} = add!(c.kinks, first(p), last(p))
function add!(ck::Kinks, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital}
    for k in ps
        add!(ck, first(k), last(k))
    end
    return ck
end
add!(c::Configuration{T}, ck::Kinks) where {T <: Orbital} = add!(c.kinks, ck...)
function add!(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital}
    add!(c1, c2.occupations)
    add!(c1, c2.kinks)
end





using CPIMC.PlaneWaves

S = sphere_with_same_spin(PlaneWave((0,0,0)),dk=1)
a = PlaneWave((-2,0,0))
b = PlaneWave((3,0,0))
c = PlaneWave((0,0,0))
d = PlaneWave((1,0,0))
e = PlaneWave((5,9,9))
f = PlaneWave((1,1,1))

g = PlaneWave(a.vec + e.vec - f.vec, Up)
h = PlaneWave(c.vec + d.vec - g.vec, Up)

i = PlaneWave((39,100,-5))
j = PlaneWave((100,0,6))

sd = Kinks( ImgTime(0.2) => T4(a,b,c,d),
       ImgTime(0.5) => T4(c,d,a,b),
       ImgTime(0.6) => T4(b,a,d,c),
       ImgTime(0.8) => T4(d,c,b,a) )

conf = Configuration(sphere(PlaneWave((0,0,0),Up),dk=4),sd)

droporbs = (c,d,f)
addorbs = (e,i,j)
dropkinks = (ImgTime(0.2) => T4(a,b,c,d), ImgTime(0.5) => T4(c,d,a,b))
addkinks = (ImgTime(0.9) => T2(a,b), ImgTime(0.1) => T4(a,a,a,a))

droporbs_set = Set((c,d,f))
addorbs_set = Set((e,i,j))
dropkinks_set = Set( (ImgTime(0.2) => T4(a,b,c,d), ImgTime(0.5) => T4(c,d,a,b)) )
addkinks_set = Set( (ImgTime(0.9) => T2(a,b), ImgTime(0.1) => T4(a,a,a,a)) )

using BenchmarkTools

println("\t Benchmark ")
println(repeat("=", 80))

println("\t construct empty Step:")
println(repeat(" -", 40), "\n")

print(typeof(StepA()), "\t:\t")
@btime StepA();
println("\n")
print(typeof(StepB()), "\t:\t")
@btime StepB();
println("\n")
print(typeof(StepC()), "\t:\t")
@btime StepC();
println("\n")
print(typeof(emptyStepD(Nothing)), "\t:\t")
@btime emptyStepD(Nothing);
println("\n")
print(typeof(emptyStepE(Nothing)), "\t:\t")
@btime emptyStepE(Nothing);
println("\n")
print(typeof(StepF()), "\t:\t")
@btime StepF();
println("\n")
print(typeof(StepG()), "\t:\t")
@btime StepG();
println("\n")
print("old version \t:\t")
@btime Step();
println("\n")


println("\t construct nonempty Step:\t")
print("drop ", length(droporbs), " orbs and ", length(dropkinks), " kinks \t")
println("add ", length(addorbs), " orbs and ", length(addkinks), " kinks")
println(repeat(" -", 40), "\n")
println("\n")

print(StepA, "\t:\t")
@btime StepA(droporbs,dropkinks,addorbs,addkinks);
println("\n")
print(StepB, "\t:\t")
@btime StepB(droporbs,dropkinks,addorbs,addkinks);
println("\n")
print(StepC, "\t:\t")
@btime StepC(droporbs,dropkinks,addorbs,addkinks);
println("\n")
print(StepD, "\t:\t")
@btime StepD(droporbs,dropkinks,addorbs,addkinks);
println("\n")
print(StepE, "\t:\t")
@btime StepE(droporbs,dropkinks,addorbs,addkinks);
println("\n")
print(StepF, "\t:\t")
@btime StepF(droporbs,dropkinks,addorbs,addkinks);
println("\n")
print(StepG, "\t:\t")
@btime StepG(droporbs,dropkinks,addorbs,addkinks);
println("\n")
print("old version \t:\t")
@btime Step(Configuration(Set(droporbs), dropkinks...), Configuration(Set(addorbs), addkinks...));
println("\n")


print("\t apply_step!:\t")
print("drop ", length(droporbs), " orbs and ", length(dropkinks), " kinks \t")
println("add ", length(addorbs), " orbs and ", length(addkinks), " kinks")
println(repeat(" -", 40), "\n")

x = StepA(droporbs,dropkinks,addorbs,addkinks);
print(typeof(x), "\t:\t")
@btime apply_step!(co,$x) setup=(co=$conf);
println("\n")
x = StepB(droporbs,dropkinks,addorbs,addkinks);
print(typeof(x), "\t:\t")
@btime apply_step!(co,$x) setup=(co=$conf);
println("\n")
x = StepC(droporbs,dropkinks,addorbs,addkinks);
print(typeof(x), "\t:\t")
@btime apply_step!(co,$x) setup=(co=$conf);
println("\n")
x = StepD(droporbs,dropkinks,addorbs,addkinks);
print(typeof(x), "\t:\t")
@btime apply_step!(co,$x) setup=(co=$conf);
println("\n")
x = StepE(droporbs,dropkinks,addorbs,addkinks);
print(typeof(x), "\t:\t")
@btime apply_step!(co,$x) setup=(co=$conf);
println("\n")
x = StepF(droporbs,dropkinks,addorbs,addkinks);
print(typeof(x), "\t:\t")
@btime apply_step!(co,$x) setup=(co=$conf);
println("\n")
x = StepG(droporbs,dropkinks,addorbs,addkinks);
print(typeof(x), "\t:\t")
@btime apply_step!(co,$x) setup=(co=$conf);
println("\n")
x = Step(Configuration(Set(droporbs), dropkinks...), Configuration(Set(addorbs), addkinks...));
print("old version \t:\t")
@btime apply_step!(co,$x) setup=(co=$conf);
println("\n")
