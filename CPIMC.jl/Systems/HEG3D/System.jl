## Calculate a mapping from orbital indices to orbital quantum numbers
## e.g. kvectors
## The units are so that the possible values of each k-vector component are integer numbers (i.e. divide by 2pi/boxarea compared to a.u.)

using StaticArrays
using DataFrames
import LinearAlgebra: dot

" single particle energy for momentum vector k "
function get_energy(o::Orbital)
    dot(o.vec,o.vec)
end


function get_orbshell(o::Orbital;dw::Int=2)
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


function get_sphere(o::Orbital; dk::Int=2) :: Set{Orbital{3}}
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


### Estimators
function Ekin(c::Configuration)
    sum(get_energy(n) for n in c.occupations)
end

function occVec(c::Configuration)
    get_energy.(c.occupations)
end

function particleNumber(c::Configuration)
  return length(c.occupations)
end

### Units
function get_beta_internal(theta, N)
  return ((2*pi)^2)/(((6*(pi^2)*N)^(2/3))*theta)
end
