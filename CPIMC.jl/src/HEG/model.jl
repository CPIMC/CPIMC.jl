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
struct Orbital{D}
    "D-dimensional excitation vector"
    vec :: StaticVector{D,Int}
    "spin"
    spin :: Int
end

Orbital(v::Tuple,s=0) = Orbital(SVector(v),s)


function get_beta_internal(theta, N)
  return ((2*pi)^2)/(((6*(pi^2)*N)^(2/3))*theta)
end

" single particle energy for momentum vector k "
function get_energy(o::Orbital)
    dot(o.vec,o.vec)
end


function get_orbshell(o::Orbital{1};dw::Int=2)
    eq = get_energy(o)
    qmax = Int(floor(eq))

    os = Set{Orbital{D}}()

    for x in -qmax:qmax
        if abs(x*x - eq) <= dw
            push!(os, Orbital((x)))
        end
    end

    os
end

function get_orbshell(o::Orbital{1};dw::Int=2)
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

function get_sphere(o::Orbital{1}; dk::Int=2)
    os = Set{Orbital{1}}()

    for x in -dk:dk
        if x*x <= dk*dk
            push!(os, Orbital(o.vec+SVector(x),0))
        end
    end

    os
end

function get_sphere(o::Orbital{2}; dk::Int=2)
    os = Set{Orbital{2}}()

    for x in -dk:dk
        for y in -dk:dk
            if x*x + y*y <= dk*dk
                push!(os, Orbital(o.vec+SVector(x,y),0))
            end
        end
    end

    os
end

function get_sphere(o::Orbital{3}; dk::Int=2)
    os = Set{Orbital{3}}()

    for x in -dk:dk
        for y in -dk:dk
            for z in -dk:dk
                if x*x + y*y + z*z <= dk*dk
                    push!(os, Orbital(o.vec+SVector(x,y,z),0))
                end
            end
        end
    end

    os
end
