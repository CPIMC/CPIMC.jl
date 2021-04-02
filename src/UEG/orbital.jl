using StaticArrays

import LinearAlgebra: dot

@enum Spin Down Up

"""
    flip(s::Spin)

return opposite spin projection
"""
function flip(s::Spin)
    if s == Down
        Up
    else
        Down
    end
end


"representation of single particle state"
struct OrbitalHEG{D} <: Orbital
    "D-dimensional momentum vector"
    vec :: StaticVector{D,Int}
    "spin"
    spin :: Spin
end

OrbitalHEG(v::Tuple,s=Up) = OrbitalHEG(SVector(v),s)

" return the energy of a single diagonal single-particle matrix element "
function energy(o::OrbitalHEG)
    dot(o.vec,o.vec)
end

"""
    w(i::Orbital, j::Orbital, k::Orbital, l::Orbital)

return the two-particle matrix element of the Coulomb interaction
zero is returned if momentum or spin is not conserved, in consequence of the Bloch-theorem and the spin kronecker-delta in the plane-spin-wave basis

the singular contribution (i.vec == k.vec) is also removed, i.e. set to zero
this property is not a result of the Bloch theorem for the two-particle Coulomb matrix element in the plane wave basis
but is added here to allow for straightforward summation over all combination of orbitals without explictly excluding this component
"""
function w(i::Orbital, j::Orbital, k::Orbital, l::Orbital)
    if !iszero(i.vec + j.vec - k.vec - l.vec) | (i.spin != k.spin) | (j.spin != l.spin)# momentum and spin conservation
        return 0.0
    elseif i.vec == k.vec
        # here return 0 in order to remove this term from sums over all occupations since it is canceled by the uniform background in the TD-Limes
        0.0
    else
        return 1.0 / dot(i.vec - k.vec, i.vec - k.vec)
    end
end

"""
    flip(o::OrbitalHEG)
return an OrbitalHEG with the same momentum vector but opposite spin projection as the given o::OrbitalHEG """
flip(o::OrbitalHEG) = OrbitalHEG(tuple(o.vec...),flip(o.spin))


"""
    dimension(o::OrbitalHEG{D}) where {D}
return the dimension of the momentum vector of an Orbital"""
dimension(o::OrbitalHEG{D}) where {D} = D

"""
    dimension(os::Set{<:OrbitalHEG{D}}) where {D}
return the dimension of the momentum vectors of a Set{<:Orbital{D}}
"""
dimension(os::Set{<:OrbitalHEG{D}}) where {D} = D

function shell(o::OrbitalHEG{1};dw::Int=2)
    eq = energy(o)
    qmax = Int(floor(eq))

    os = Set{OrbitalHEG{1}}()

    for x in -qmax:qmax
        if abs(x*x - eq) <= dw
            push!(os, OrbitalHEG((x,1)))
            push!(os, OrbitalHEG((x,-1)))
        end
    end

    os
end

function shell(o::OrbitalHEG{2};dw::Int=2)
    eq = energy(o)
    qmax = Int(floor(eq))

    os = Set{OrbitalHEG{2}}()

    for x in -qmax:qmax
        for y in -qmax:qmax
            if abs(x*x + y*y - eq) <= dw
                push!(os, OrbitalHEG((x,y),1))
                push!(os, OrbitalHEG((x,y),-1))
            end
        end
    end

    os
end

function shell(o::OrbitalHEG{3};dw::Int=2)
    eq = energy(o)
    qmax = Int(floor(eq))

    os = Set{OrbitalHEG{3}}()

    for x in -qmax:qmax
        for y in -qmax:qmax
            for z in -qmax:qmax
                if abs( x*x + y*y + z*z - eq ) <= dw
                    push!(os, OrbitalHEG((x,y,z),1))
                    push!(os, OrbitalHEG((x,y,z),-1))
                end
            end
        end
    end

    os
end

function sphere_with_same_spin(o::OrbitalHEG{1}; dk::Number=2)
    os = Set{OrbitalHEG{1}}()

    for x in -Int(ceil(dk)):Int(ceil(dk))
        if x*x <= dk*dk
            push!(os, OrbitalHEG(o.vec+SVector(x),o.spin))
        end
    end
    os
end

function sphere_with_same_spin(o::OrbitalHEG{2}; dk::Number=2)
    os = Set{OrbitalHEG{2}}()

    for x in -Int(ceil(dk)):Int(ceil(dk))
        for y in -Int(ceil(dk)):Int(ceil(dk))
            if x*x + y*y <= dk*dk
                push!(os, OrbitalHEG(o.vec+SVector(x,y),o.spin))
            end
        end
    end
    os
end


function sphere_with_same_spin(o::OrbitalHEG{3}; dk::Number=2)
    os = Set{OrbitalHEG{3}}()

    for x in -Int(ceil(dk)):Int(ceil(dk))
        for y in -Int(ceil(dk)):Int(ceil(dk))
            for z in -Int(ceil(dk)):Int(ceil(dk))
                if x*x + y*y + z*z <= dk*dk
                    push!(os, OrbitalHEG(o.vec+SVector(x,y,z),o.spin))
                end
            end
        end
    end
    os
end

" return two spheres of each spin with radius dk around the wavevector of the argument orbital o"
sphere(o::OrbitalHEG; dk=2) = union(sphere_with_same_spin(o,dk=dk),sphere_with_same_spin(OrbitalHEG(o.vec,flip(o.spin)),dk=dk))

function orbs_with_spin(orbitals::Set{OrbitalHEG{D}},spin::Int) where {D}
    orbs_s = Set{OrbitalHEG{D}}()
    for orb in orbitals
        if orb.spin==spin
            push!(orbs_s, orb)
        end
    end
    return orbs_s
end
