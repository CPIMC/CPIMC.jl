using StaticArrays

import LinearAlgebra: dot

struct Ensemble
  "Brueckner parameter"
  rs :: Float64

  "reduced temperature"
  beta :: Float64

  "particle number"
  N :: Int
end

"representation of single particle state"
struct Orbital{D} <: Basis
    "D-dimensional excitation vector"
    vec :: StaticVector{D,Int}
    "spin"
    spin :: Int
end

Orbital(v::Tuple,s=0) = Orbital(SVector(v),s)

function get_beta_internal(theta, N)
  return ((2*pi)^2)/(((6*(pi^2)*N)^(2/3))*theta)
end

""" return the energy of a single diagonal single-particle matrix element
    for a fully occupied orbital """
function get_energy(o::Orbital)
    dot(o.vec,o.vec)
end

## method declarations for the off-diagonal energy given by kinks in the UEG
# TODO: It may not be good style to use the same function name for the off-diagnoal energy as for the diagonal energy unless we dispatch at some point.
#       On the other hand, this notation is semantically convenient since it is short and conceptually unambiguous (i.e. it is clear what is the difference between the energy of an orbital (diagonal) and the energy of a kink (off-diagonal)).

""" return the energy of a single off-diagonal single-particle matrix element (a.k.a. kink)
    for a two-particle excitation (i.e. T4)"""
function get_energy(t4::T4)
    @assert iszero(t4.i + t4.j + t4.k + t4.l) "Momentum conservation violated by kink $(t4[2]). All excitations must conserve the total momentum in the UEG."
    @assert !iszero(t4.i - t4.k) "Divergence in off-diagonal single-particle matrix element. All excitations must conserve the total momentum in the UEG."
    1.0 / abs( t4.i - t4.k)
end

""" return the energy of a single off-diagonal matrix element (a.k.a. kink)
    for a one-particle excitation (i.e. T2)"""
function get_energy(t2::T2)
    throw(ErrorException("There are no T2-kinks in the UEG."))
end


""" return the energy of a configuration """
function get_energy(c::Configuration)
end

function get_orbshell(o::Orbital{1};dw::Int=2)
    eq = get_energy(o)
    qmax = Int(floor(eq))

    os = Set{Orbital{1}}()

    for x in -qmax:qmax
        if abs(x*x - eq) <= dw
            push!(os, Orbital((x)))
        end
    end

    os
end

function get_orbshell(o::Orbital{2};dw::Int=2)
    eq = get_energy(o)
    qmax = Int(floor(eq))

    os = Set{Orbital{2}}()

    for x in -qmax:qmax
        for y in -qmax:qmax
            if abs(x*x + y*y - eq) <= dw
                push!(os, Orbital((x,y)))
            end
        end
    end

    os
end

function get_orbshell(o::Orbital{3};dw::Int=2)
    eq = get_energy(o)
    qmax = Int(floor(eq))

    os = Set{Orbital{3}}()

    for x in -qmax:qmax
        for y in -qmax:qmax
            for z in -qmax:qmax
                if abs( x*x + y*y + z*z - eq ) <= dw
                    push!(os, Orbital((x,y,z)))
                end
            end
        end
    end

    os
end

function get_sphere(o::Orbital{1}; dε::Int=4)
    os = Set{Orbital{1}}()

    for x in -dε:dε
        if x*x <= dε
            push!(os, Orbital(o.vec+SVector(x),0))
        end
    end

    os
end

function get_sphere(o::Orbital{2}; dε::Int=4)
    os = Set{Orbital{2}}()

    for x in -dε:dε
        for y in -dε:dε
            if x*x + y*y <= dε
                push!(os, Orbital(o.vec+SVector(x,y),0))
            end
        end
    end

    os
end

function get_sphere(o::Orbital{3}; dε::Int=4)
    os = Set{Orbital{3}}()

    for x in -dε:dε
        for y in -dε:dε
            for z in -dε:dε
                if x*x + y*y + z*z <= dε
                    push!(os, Orbital(o.vec+SVector(x,y,z),0))
                end
            end
        end
    end

    os
end
