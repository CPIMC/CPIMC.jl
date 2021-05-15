"""
Plane wave single-particle basis states.
"""
module PlaneWaves

export Spin, Down, Up, PlaneWave, dimension, flip, sphere, sphere_with_same_spin, ξ, fractional_spin_polarization, find_fourth_orb_for_kink

using ..CPIMC
using StaticArrays

"""
Type for storing spin projection of particles with spin 1/2.
"""
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


"""
    struct PlaneWave{D} <: Orbital

Representation of plane wave single-particle states with dimensionality `D`.

**Fields**
- `vec  :: StaticVector{D,Int}` -- momentum vector
- `spin :: Spin`                -- spin

"""
struct PlaneWave{D} <: Orbital
    vec :: SVector{D,Int}
    spin :: Spin
end

"""
    PlaneWave(v::Tuple,s=Up)

Outer constructor with default value `spin = Up`.
"""
PlaneWave(v::Tuple,s=Up) = PlaneWave(SVector(v),s)


"""
    flip(o::PlaneWave)
return an PlaneWave with the same momentum vector but opposite spin projection as the given o::PlaneWave """
flip(o::PlaneWave) = PlaneWave(tuple(o.vec...),flip(o.spin))


"""
    dimension(o::PlaneWave{D}) where {D}
return the dimension of the momentum vector of an Orbital"""
dimension(o::PlaneWave{D}) where {D} = D

"""
    dimension(os::Set{<:PlaneWave{D}}) where {D}
return the dimension of the momentum vectors of a Set{<:Orbital{D}}
"""
dimension(os::Set{<:PlaneWave{D}}) where {D} = D

"""
    sphere_with_same_spin(o::PlaneWave{1}; dk=2)
    sphere_with_same_spin(o::PlaneWave{2}; dk=2)
    sphere_with_same_spin(o::PlaneWave{3}; dk=2)

Return the set of all PlaneWaves in a sphere of radius `dk` around `o` in momentum space that also have the same spin as `o`.
"""
function sphere_with_same_spin(o::PlaneWave{1}; dk::Number=2)
    os = Set{PlaneWave{1}}()

    for x in -Int(ceil(dk)):Int(ceil(dk))
        if x*x <= dk*dk
            push!(os, PlaneWave(o.vec+SVector(x),o.spin))
        end
    end
    os
end

function sphere_with_same_spin(o::PlaneWave{2}; dk::Number=2)
    os = Set{PlaneWave{2}}()

    for x in -Int(ceil(dk)):Int(ceil(dk))
        for y in -Int(ceil(dk)):Int(ceil(dk))
            if x*x + y*y <= dk*dk
                push!(os, PlaneWave(o.vec+SVector(x,y),o.spin))
            end
        end
    end
    os
end


function sphere_with_same_spin(o::PlaneWave{3}; dk::Number=2)
    os = Set{PlaneWave{3}}()

    for x in -Int(ceil(dk)):Int(ceil(dk))
        for y in -Int(ceil(dk)):Int(ceil(dk))
            for z in -Int(ceil(dk)):Int(ceil(dk))
                if x*x + y*y + z*z <= dk*dk
                    push!(os, PlaneWave(o.vec+SVector(x,y,z),o.spin))
                end
            end
        end
    end
    os
end

"""
    sphere(o::PlaneWave{1}; dk=2)
    sphere(o::PlaneWave{2}; dk=2)
    sphere(o::PlaneWave{3}; dk=2)

Return the set of all PlaneWaves in a sphere of radius `dk` around `o` in momentum space.
"""
sphere(o::PlaneWave; dk=2) = union(sphere_with_same_spin(o,dk=dk),sphere_with_same_spin(PlaneWave(o.vec,flip(o.spin)),dk=dk))


"""
    fractional_spin_polarization(occ::Set{PlaneWave)
calculate the fractional spin polarization from a set of Orbitals which each must have a field `spin` of type `@enum Spin Down Up`,
the fractional spin polarization is defined as the ratio of the absolute difference of the number of particles which Spin Down and Spin Up and the particle number,
thus here no convention is made as to which spin component should be occupied more often, as opposed to some literature where N↑ > N↓ is used
"""
function fractional_spin_polarization(occ::Set{<:PlaneWave})
    N_up = length(filter(x -> x.spin == Up, occ))
    N_down = length(filter(x -> x.spin == Down, occ))
    N = length(occ)
    return abs(N_up - N_down) / N
end

"""
    ξ(occ::Set{PlaneWave})
calculate the fractional spin polarization from a set of Orbitals which each must have a field `spin` of type `@enum Spin Down Up`,
the fractional spin polarization is defined as the ratio of the absolute difference of the number of particles which Spin Down and Spin Up and the particle number,
thus here no convention is made as to which spin component should be occupied more often, as opposed to some literature where N↑ > N↓ is used
"""
ξ(occ::Set{<:PlaneWave}) = fractional_spin_polarization(occ)



"""
    fractional_spin_polarization(c::Configuration{<:PlaneWave})
calculate the fractional spin polarization from a set of Orbitals which each must have a field `spin` of type `@enum Spin Down Up`,
the fractional spin polarization is defined as the ratio of the absolute difference of the number of particles which Spin Down and Spin Up and the particle number,
thus here no convention is made as to which spin component should be occupied more often, as opposed to some literature where N↑ > N↓ is used
"""
function fractional_spin_polarization(c::Configuration{<:PlaneWave})# TODO: added {<:PlaneWave}
    fractional_spin_polarization(c.occupations)
end


"""
    ξ(c::Configuration{<:PlaneWave})
calculate the fractional spin polarization from a set of Orbitals which each must have a field `spin` of type `@enum Spin Down Up`,
the fractional spin polarization is defined as the ratio of the absolute difference of the number of particles which Spin Down and Spin Up and the particle number,
thus here no convention is made as to which spin component should be occupied more often, as opposed to some literature where N↑ > N↓ is used
"""
ξ(c::Configuration{<:PlaneWave}) = fractional_spin_polarization(c)


"""
    find_fourth_orb_for_kink(same_kind_ladder_operator, other_kind_ladder_operator1, other_kind_ladder_operator2)
    Find the fourth orb to build a Kink accoring to spin/momentum conservation, where the Operator thats acts on the
    fourth orb is of the same kind (creator/ annihilator) as the one that acts on same_kind_ladder_operator,
    While on other_kind_ladder_operator1 and other_kind_ladder_operator2 the other kind of operators will act.
"""
function find_fourth_orb_for_kink(same_kind_ladder_operator, other_kind_ladder_operator1, other_kind_ladder_operator2)
    @assert(in(same_kind_ladder_operator.spin, [other_kind_ladder_operator1.spin, other_kind_ladder_operator2.spin]))
    if other_kind_ladder_operator1.spin == other_kind_ladder_operator2.spin
        fourth_orb = PlaneWave(other_kind_ladder_operator1.vec + other_kind_ladder_operator2.vec - same_kind_ladder_operator.vec, same_kind_ladder_operator.spin)

    else
        fourth_orb = PlaneWave(other_kind_ladder_operator1.vec + other_kind_ladder_operator2.vec - same_kind_ladder_operator.vec, flip(same_kind_ladder_operator.spin))
    end
    return(fourth_orb)
end

end
