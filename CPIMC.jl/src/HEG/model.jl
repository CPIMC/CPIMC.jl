using StaticArrays

import LinearAlgebra: dot

struct Ensemble
  "Brueckner parameter"
  rs :: Float64
  "reduced temperature"
  β :: Float64
  "particle number"
  N :: Int
end

"representation of single particle state"
struct OrbitalHEG{D} <: Orbital
    "D-dimensional excitation vector"
    vec :: StaticVector{D,Int}
    "spin"
    spin :: Int
end

OrbitalHEG(v::Tuple,s=1) = OrbitalHEG(SVector(v),s)


function get_β_internal(θ::Float64, N::Int)
  return ((2*pi)^2)/(((6*(pi^2)*N)^(2/3))*θ)   #if using unpolarized system need a factor 1/2 under the ^(2/3)
end

function lambda(N::Int, rs::Float64)
    #Warum Faktor 2 am ende?
    return (4/((2*pi)^3)) * (((4*pi)/3)^(1/3)) * rs * N^(1/3) * 2
end

" single particle energy for momentum vector k "
function get_energy(o::OrbitalHEG)
    dot(o.vec,o.vec)
end

function E_Ry(E_internal::Float64, lam::Float64)
    return (E_internal * 16/((2*pi)^4 * (lam/2)^2))   #if we change the factor 2 in lambda we have to change the factor lam/2 in this formula
end

function abs_E_mad(N::Int, lam::Float64) #internal units
    return 2.83729747948527 * pi/2.0 * N * (lam/2)   #if we change the factor 2 in lambda we have to change the factor lam/2 in this formula
end


function get_abs_offdiagonal_element(e::Ensemble,c::Configuration,kink::T4{OrbitalHEG{3}})
    wijkl = 0
    if kink.i.spin == kink.j.spin
        wijkl +=  1/dot((kink.i.vec-kink.k.vec), (kink.i.vec-kink.k.vec)) -
                    1/dot((kink.i.vec-kink.l.vec), (kink.i.vec-kink.l.vec))
    elseif kink.i.spin == kink.k.spin
        wijkl +=  1/dot((kink.i.vec-kink.k.vec), (kink.i.vec-kink.k.vec))
    elseif kink.i.spin == kink.l.spin
        wijkl -= 1/dot((kink.i.vec-kink.l.vec), (kink.i.vec-kink.l.vec))
    else
        @assert false
    end
    #the factor lamda/2 is due to the use of internal units
    wijkl *= lambda(e.N,e.rs)/2
    # We sample with the weight of antisymmetriesed matrixelement but we do not restrict
    # the orders of indizies of our possible kinks. We therefor need an extrra factor 1/4 in the weight-funktion
    return abs(wijkl) * 1/4
end



function get_orbshell(o::OrbitalHEG{1};dw::Int=2)
    eq = get_energy(o)
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

function get_orbshell(o::OrbitalHEG{2};dw::Int=2)
    eq = get_energy(o)
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

function get_orbshell(o::OrbitalHEG{3};dw::Int=2)
    eq = get_energy(o)
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

function get_sphere_with_same_spin(o::OrbitalHEG{1}; dk::Number=2)
    os = Set{OrbitalHEG{1}}()

    for x in -Int(ceil(dk)):Int(ceil(dk))
        if x*x <= dk*dk
            push!(os, OrbitalHEG(o.vec+SVector(x),o.spin))
        end
    end
    os
end

function get_sphere_with_same_spin(o::OrbitalHEG{2}; dk::Number=2)
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


function get_sphere_with_same_spin(o::OrbitalHEG{3}; dk::Number=2)
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

function get_orbs_with_spin(orbitals::Set{OrbitalHEG{3}},spin::Int)
    orbs_s = Set{OrbitalHEG{3}}()
    for orb in orbitals
        if orb.spin==spin
            push!(orbs_s, orb)
        end
    end
    return orbs_s
end



#calculates the change in the diagonal interaction when changing ocupation between τ1 and τ2 accoring to left_kink
#already multiplied by e.β
function get_change_diagonal_interaction(c::Configuration, e::Ensemble, left_kink::T4, τ1, τ2)
    orb_a = left_kink.i
    orb_b = left_kink.j
    orb_c = left_kink.k
    orb_d = left_kink.l
    delta_τ12 = τ2 - τ1
    if delta_τ12 < 0
        delta_τ12 += 1
    end
    delta_di = delta_τ12 * (lambda(e.N,e.rs)/2) * (1/dot((orb_a.vec-orb_b.vec),(orb_a.vec-orb_b.vec)) -
                                        1/dot((orb_c.vec-orb_d.vec),(orb_c.vec-orb_d.vec)))
    occs = get_occupations_at(c, τ1)
    for occ in occs
        if ((occ == orb_c) | (occ == orb_d) | (occ == orb_a) | (occ == orb_b))
            "nix"
        else
            delta_di += delta_τ12 * (lambda(e.N,e.rs)/2) * (1/dot((occ.vec-orb_a.vec),(occ.vec-orb_a.vec))
                                                + 1/dot((occ.vec-orb_b.vec),(occ.vec-orb_b.vec))
                                                - 1/dot((occ.vec-orb_c.vec),(occ.vec-orb_c.vec))
                                                - 1/dot((occ.vec-orb_d.vec),(occ.vec-orb_d.vec)))
        end
    end
    if length(c.kinks) == 0
        return -delta_di * e.β
    end
    kink_semi_token = searchsortedfirst(c.kinks,τ1)
    if kink_semi_token == pastendsemitoken(c.kinks)
        kink_semi_token = startof(c.kinks)
    end
    τ_kink,kink = deref((c.kinks,kink_semi_token))
    #The kink at τ1 is already considered in occs
    if τ_kink == τ1
        kink_semi_token = advance((c.kinks,kink_semi_token))
        if kink_semi_token == pastendsemitoken(c.kinks)
            kink_semi_token = startof(c.kinks)
        end
        τ_kink,kink = deref((c.kinks,kink_semi_token))
    end
    loop_counter = 0
    while ((τ1 < τ_kink < τ2) | (τ_kink < τ2 < τ1) | (τ2 < τ1 < τ_kink)) & (loop_counter < length(c.kinks))
        delta_τ = τ2 - τ_kink
        if delta_τ < 0
            delta_τ += 1
        end

        delta_di += delta_τ * (lambda(e.N,e.rs)/2) *
                                (1/dot((kink.i.vec-orb_a.vec),(kink.i.vec-orb_a.vec))
                                    + 1/dot((kink.i.vec-orb_b.vec),(kink.i.vec-orb_b.vec))
                                    - 1/dot((kink.i.vec-orb_c.vec),(kink.i.vec-orb_c.vec))
                                    - 1/dot((kink.i.vec-orb_d.vec),(kink.i.vec-orb_d.vec)))
        delta_di += delta_τ * (lambda(e.N,e.rs)/2) *
                                (1/dot((kink.j.vec-orb_a.vec),(kink.j.vec-orb_a.vec))
                                    + 1/dot((kink.j.vec-orb_b.vec),(kink.j.vec-orb_b.vec))
                                    - 1/dot((kink.j.vec-orb_c.vec),(kink.j.vec-orb_c.vec))
                                    - 1/dot((kink.j.vec-orb_d.vec),(kink.j.vec-orb_d.vec)))
        delta_di -= delta_τ * (lambda(e.N,e.rs)/2) *
                                (1/dot((kink.k.vec-orb_a.vec),(kink.k.vec-orb_a.vec))
                                    + 1/dot((kink.k.vec-orb_b.vec),(kink.k.vec-orb_b.vec))
                                    - 1/dot((kink.k.vec-orb_c.vec),(kink.k.vec-orb_c.vec))
                                    - 1/dot((kink.k.vec-orb_d.vec),(kink.k.vec-orb_d.vec)))
        delta_di -= delta_τ * (lambda(e.N,e.rs)/2) *
                                (1/dot((kink.l.vec-orb_a.vec),(kink.l.vec-orb_a.vec))
                                    + 1/dot((kink.l.vec-orb_b.vec),(kink.l.vec-orb_b.vec))
                                    - 1/dot((kink.l.vec-orb_c.vec),(kink.l.vec-orb_c.vec))
                                    - 1/dot((kink.l.vec-orb_d.vec),(kink.l.vec-orb_d.vec)))

        kink_semi_token = advance((c.kinks,kink_semi_token))
        if kink_semi_token == pastendsemitoken(c.kinks)
            kink_semi_token = startof(c.kinks)
        end
        τ_kink,kink = deref((c.kinks,kink_semi_token))
        loop_counter += 1
    end
    return -delta_di * e.β
end

#This function assumes for now that all kinks in configuration.kinks are of type 4
function get_change_diagonal_interaction(c::Configuration, e::Ensemble, left_kink::T2, τ1, τ2)
    orb_a = left_kink.i
    orb_b = left_kink.j
    delta_τ12 = τ2 - τ1
    if delta_τ12 < 0
        delta_τ12 += 1
    end
    delta_di = 0
    occs = get_occupations_at(c, τ1)
    @assert(!in(orb_a, occs))
    for occ in occs
        if ((occ == orb_b) | (occ == orb_a))
            "nix"
        else
            delta_di += delta_τ12 * (lambda(e.N,e.rs)/2) * (1/dot((occ.vec-orb_a.vec),(occ.vec-orb_a.vec))
                                                - 1/dot((occ.vec-orb_b.vec),(occ.vec-orb_b.vec)))

        end
    end
    if length(c.kinks) == 0
        return -delta_di * e.β
    end
    kink_semi_token = searchsortedfirst(c.kinks,τ1)
    if kink_semi_token == pastendsemitoken(c.kinks)
        kink_semi_token = startof(c.kinks)
    end
    τ_kink,kink = deref((c.kinks,kink_semi_token))
    #The kink at τ1 is already considered in occs
    if τ_kink == τ1
        kink_semi_token = advance((c.kinks,kink_semi_token))
        if kink_semi_token == pastendsemitoken(c.kinks)
            kink_semi_token = startof(c.kinks)
        end
        τ_kink,kink = deref((c.kinks,kink_semi_token))
    end
    loop_counter = 0
    while ((τ1 < τ_kink < τ2) | (τ_kink < τ2 < τ1) | (τ2 < τ1 < τ_kink)) & (loop_counter < length(c.kinks))
        delta_τ = τ2 - τ_kink
        if delta_τ < 0
            delta_τ += 1
        end

        delta_di += delta_τ * (lambda(e.N,e.rs)/2) *
                                (1/dot((kink.i.vec-orb_a.vec),(kink.i.vec-orb_a.vec))
                                    - 1/dot((kink.i.vec-orb_b.vec),(kink.i.vec-orb_b.vec)))
        delta_di += delta_τ * (lambda(e.N,e.rs)/2) *
                                (1/dot((kink.j.vec-orb_a.vec),(kink.j.vec-orb_a.vec))
                                    - 1/dot((kink.j.vec-orb_b.vec),(kink.j.vec-orb_b.vec)))
        delta_di -= delta_τ * (lambda(e.N,e.rs)/2) *
                                (1/dot((kink.k.vec-orb_a.vec),(kink.k.vec-orb_a.vec))
                                    - 1/dot((kink.k.vec-orb_b.vec),(kink.k.vec-orb_b.vec)))
        delta_di -= delta_τ * (lambda(e.N,e.rs)/2) *
                                (1/dot((kink.l.vec-orb_a.vec),(kink.l.vec-orb_a.vec))
                                    - 1/dot((kink.l.vec-orb_b.vec),(kink.l.vec-orb_b.vec)))
        kink_semi_token = advance((c.kinks,kink_semi_token))
        if kink_semi_token == pastendsemitoken(c.kinks)
            kink_semi_token = startof(c.kinks)
        end
        τ_kink,kink = deref((c.kinks,kink_semi_token))
        loop_counter += 1
    end
    return -delta_di * e.β
end
