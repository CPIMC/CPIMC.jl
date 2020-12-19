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
struct Orbital_HEG{D} <: Orbital
    "D-dimensional excitation vector"
    vec :: StaticVector{D,Int}
    "spin"
    spin :: Int
end

Orbital_HEG(v::Tuple,s=1) = Orbital_HEG(SVector(v),s)


function get_beta_internal(theta::Float64, N::Int)
  return ((2*pi)^2)/(((6*(pi^2)*N)^(2/3))*theta)   #if using unpolarized system need a factor 1/2 under the ^(2/3)
end

function lambda(N::Int, rs::Float64)
    #Warum Faktor 2 am ende?
    return (4/((2*pi)^3)) * (((4*pi)/3)^(1/3)) * rs * N^(1/3) * 2
end

" single particle energy for momentum vector k "
function get_energy(o::Orbital_HEG)
    dot(o.vec,o.vec)
end

function E_Ry(E_internal::Float64, lam::Float64)
    return (E_internal * 16/((2*pi)^4 * (lam/2)^2))   #if we change the factor 2 in lambda we have to change the factor lam/2 in this formula
end

function abs_E_Mad(N::Int, lam::Float64) #internal units
    return 2.83729747948527 * pi/2.0 * N * (lam/2)   #if we change the factor 2 in lambda we have to change the factor lam/2 in this formula
end


function get_abs_offdiagonal_element(e::Ensemble,c::Configuration,Kink::T4{Orbital_HEG{3}})
    wijkl = 0
    if Kink.i.spin == Kink.j.spin
        wijkl +=  1/dot((Kink.i.vec-Kink.k.vec), (Kink.i.vec-Kink.k.vec)) -
                    1/dot((Kink.i.vec-Kink.l.vec), (Kink.i.vec-Kink.l.vec))
    elseif Kink.i.spin == Kink.k.spin
        wijkl +=  1/dot((Kink.i.vec-Kink.k.vec), (Kink.i.vec-Kink.k.vec))
    elseif Kink.i.spin == Kink.l.spin
        wijkl -= 1/dot((Kink.i.vec-Kink.l.vec), (Kink.i.vec-Kink.l.vec))
    else
        @assert false
    end
    #the factor lamda/2 is due to the use of internal units
    wijkl *= lambda(e.N,e.rs)/2
    # We sample with the weight of antisymmetriesed matrixelement but we do not restrict
    # the orders of indizies of our possible kinks. We therefor need an extrra factor 1/4 in the weight-funktion
    return abs(wijkl) * 1/4
end



function get_orbshell(o::Orbital_HEG{1};dw::Int=2)
    eq = get_energy(o)
    qmax = Int(floor(eq))

    os = Set{Orbital_HEG{1}}()

    for x in -qmax:qmax
        if abs(x*x - eq) <= dw
            push!(os, Orbital_HEG((x,1)))
            push!(os, Orbital_HEG((x,-1)))
        end
    end

    os
end

function get_orbshell(o::Orbital_HEG{2};dw::Int=2)
    eq = get_energy(o)
    qmax = Int(floor(eq))

    os = Set{Orbital_HEG{2}}()

    for x in -qmax:qmax
        for y in -qmax:qmax
            if abs(x*x + y*y - eq) <= dw
                push!(os, Orbital_HEG((x,y),1))
                push!(os, Orbital_HEG((x,y),-1))
            end
        end
    end

    os
end

function get_orbshell(o::Orbital_HEG{3};dw::Int=2)
    eq = get_energy(o)
    qmax = Int(floor(eq))

    os = Set{Orbital_HEG{3}}()

    for x in -qmax:qmax
        for y in -qmax:qmax
            for z in -qmax:qmax
                if abs( x*x + y*y + z*z - eq ) <= dw
                    push!(os, Orbital_HEG((x,y,z),1))
                    push!(os, Orbital_HEG((x,y,z),-1))
                end
            end
        end
    end

    os
end

function get_sphere(o::Orbital_HEG{1}; dk::Int=2)
    os = Set{Orbital_HEG{1}}()

    for x in -dk:dk
        if x*x <= dk*dk
            push!(os, Orbital_HEG(o.vec+SVector(x),1))
            push!(os, Orbital_HEG(o.vec+SVector(x),-1))
        end
    end
    os
end

function get_sphere(o::Orbital_HEG{2}; dk::Int=2)
    os = Set{Orbital_HEG{2}}()

    for x in -dk:dk
        for y in -dk:dk
            if x*x + y*y <= dk*dk
                push!(os, Orbital_HEG(o.vec+SVector(x,y),1))
                push!(os, Orbital_HEG(o.vec+SVector(x,y),-1))
            end
        end
    end
    os
end


function get_sphere(o::Orbital_HEG{3}; dk::Int=2)
    os = Set{Orbital_HEG{3}}()

    for x in -dk:dk
        for y in -dk:dk
            for z in -dk:dk
                if x*x + y*y + z*z <= dk*dk
                    push!(os, Orbital_HEG(o.vec+SVector(x,y,z),1))
                    push!(os, Orbital_HEG(o.vec+SVector(x,y,z),-1))
                end
            end
        end
    end
    os
end

function get_orbs_with_spin(orbitals::Set{Orbital_HEG{3}},spin::Int)
    orbs_s = Set{Orbital_HEG{3}}()
    for orb in orbitals
        if orb.spin==spin
            push!(orbs_s, orb)
        end
    end
    return orbs_s
end



#calculates the change in the diagonal interaction when changing ocupation between Tau1 and Tau2 accoring to LeftKink
#already multiplied by e.beta
function get_change_diagonal_interaction(c::Configuration, e::Ensemble, LeftKink::T4, Tau1, Tau2)
    orb_a = LeftKink.i
    orb_b = LeftKink.j
    orb_c = LeftKink.k
    orb_d = LeftKink.l
    delta_Tau12 = Tau2 - Tau1
    if delta_Tau12 < 0
        delta_Tau12 += 1
    end
    delta_di = delta_Tau12 * (lambda(e.N,e.rs)/2) * (1/dot((orb_a.vec-orb_b.vec),(orb_a.vec-orb_b.vec)) -
                                        1/dot((orb_c.vec-orb_d.vec),(orb_c.vec-orb_d.vec)))
    occs = get_occupations_at(c, Tau1)
    for occ in occs
        if ((occ == orb_c) | (occ == orb_d) | (occ == orb_a) | (occ == orb_b))
            "nix"
        else
            delta_di += delta_Tau12 * (lambda(e.N,e.rs)/2) * (1/dot((occ.vec-orb_a.vec),(occ.vec-orb_a.vec))
                                                + 1/dot((occ.vec-orb_b.vec),(occ.vec-orb_b.vec))
                                                - 1/dot((occ.vec-orb_c.vec),(occ.vec-orb_c.vec))
                                                - 1/dot((occ.vec-orb_d.vec),(occ.vec-orb_d.vec)))
        end
    end
    if length(c.kinks) == 0
        return -delta_di * e.beta
    end
    Kink_semi_token = searchsortedfirst(c.kinks,Tau1)
    if Kink_semi_token == pastendsemitoken(c.kinks)
        Kink_semi_token = startof(c.kinks)
    end
    Tau_Kink,Kink = deref((c.kinks,Kink_semi_token))
    #The Kink at Tau1 is already considered in occs
    if Tau_Kink == Tau1
        Kink_semi_token = advance((c.kinks,Kink_semi_token))
        if Kink_semi_token == pastendsemitoken(c.kinks)
            Kink_semi_token = startof(c.kinks)
        end
        Tau_Kink,Kink = deref((c.kinks,Kink_semi_token))
    end
    loop_counter = 0
    while ((Tau1 < Tau_Kink < Tau2) | (Tau_Kink < Tau2 < Tau1) | (Tau2 < Tau1 < Tau_Kink)) & (loop_counter < length(c.kinks))
        delta_Tau = Tau2 - Tau_Kink
        if delta_Tau < 0
            delta_Tau += 1
        end

        delta_di += delta_Tau * (lambda(e.N,e.rs)/2) *
                                (1/dot((Kink.i.vec-orb_a.vec),(Kink.i.vec-orb_a.vec))
                                    + 1/dot((Kink.i.vec-orb_b.vec),(Kink.i.vec-orb_b.vec))
                                    - 1/dot((Kink.i.vec-orb_c.vec),(Kink.i.vec-orb_c.vec))
                                    - 1/dot((Kink.i.vec-orb_d.vec),(Kink.i.vec-orb_d.vec)))
        delta_di += delta_Tau * (lambda(e.N,e.rs)/2) *
                                (1/dot((Kink.j.vec-orb_a.vec),(Kink.j.vec-orb_a.vec))
                                    + 1/dot((Kink.j.vec-orb_b.vec),(Kink.j.vec-orb_b.vec))
                                    - 1/dot((Kink.j.vec-orb_c.vec),(Kink.j.vec-orb_c.vec))
                                    - 1/dot((Kink.j.vec-orb_d.vec),(Kink.j.vec-orb_d.vec)))
        delta_di -= delta_Tau * (lambda(e.N,e.rs)/2) *
                                (1/dot((Kink.k.vec-orb_a.vec),(Kink.k.vec-orb_a.vec))
                                    + 1/dot((Kink.k.vec-orb_b.vec),(Kink.k.vec-orb_b.vec))
                                    - 1/dot((Kink.k.vec-orb_c.vec),(Kink.k.vec-orb_c.vec))
                                    - 1/dot((Kink.k.vec-orb_d.vec),(Kink.k.vec-orb_d.vec)))
        delta_di -= delta_Tau * (lambda(e.N,e.rs)/2) *
                                (1/dot((Kink.l.vec-orb_a.vec),(Kink.l.vec-orb_a.vec))
                                    + 1/dot((Kink.l.vec-orb_b.vec),(Kink.l.vec-orb_b.vec))
                                    - 1/dot((Kink.l.vec-orb_c.vec),(Kink.l.vec-orb_c.vec))
                                    - 1/dot((Kink.l.vec-orb_d.vec),(Kink.l.vec-orb_d.vec)))

        Kink_semi_token = advance((c.kinks,Kink_semi_token))
        if Kink_semi_token == pastendsemitoken(c.kinks)
            Kink_semi_token = startof(c.kinks)
        end
        Tau_Kink,Kink = deref((c.kinks,Kink_semi_token))
        loop_counter += 1
    end
    return -delta_di * e.beta
end

#This function assumes for now that all kinks in configuration.kinks are of type 4
function get_change_diagonal_interaction(c::Configuration, e::Ensemble, LeftKink::T2, Tau1, Tau2)
    orb_a = LeftKink.i
    orb_b = LeftKink.j
    delta_Tau12 = Tau2 - Tau1
    if delta_Tau12 < 0
        delta_Tau12 += 1
    end
    delta_di = 0
    occs = get_occupations_at(c, Tau1)
    @assert(!in(orb_a, occs))
    for occ in occs
        if ((occ == orb_b) | (occ == orb_a))
            "nix"
        else
            delta_di += delta_Tau12 * (lambda(e.N,e.rs)/2) * (1/dot((occ.vec-orb_a.vec),(occ.vec-orb_a.vec))
                                                - 1/dot((occ.vec-orb_b.vec),(occ.vec-orb_b.vec)))

        end
    end
    if length(c.kinks) == 0
        return -delta_di * e.beta
    end
    Kink_semi_token = searchsortedfirst(c.kinks,Tau1)
    if Kink_semi_token == pastendsemitoken(c.kinks)
        Kink_semi_token = startof(c.kinks)
    end
    Tau_Kink,Kink = deref((c.kinks,Kink_semi_token))
    #The Kink at Tau1 is already considered in occs
    if Tau_Kink == Tau1
        Kink_semi_token = advance((c.kinks,Kink_semi_token))
        if Kink_semi_token == pastendsemitoken(c.kinks)
            Kink_semi_token = startof(c.kinks)
        end
        Tau_Kink,Kink = deref((c.kinks,Kink_semi_token))
    end
    loop_counter = 0
    while ((Tau1 < Tau_Kink < Tau2) | (Tau_Kink < Tau2 < Tau1) | (Tau2 < Tau1 < Tau_Kink)) & (loop_counter < length(c.kinks))
        delta_Tau = Tau2 - Tau_Kink
        if delta_Tau < 0
            delta_Tau += 1
        end

        delta_di += delta_Tau * (lambda(e.N,e.rs)/2) *
                                (1/dot((Kink.i.vec-orb_a.vec),(Kink.i.vec-orb_a.vec))
                                    - 1/dot((Kink.i.vec-orb_b.vec),(Kink.i.vec-orb_b.vec)))
        delta_di += delta_Tau * (lambda(e.N,e.rs)/2) *
                                (1/dot((Kink.j.vec-orb_a.vec),(Kink.j.vec-orb_a.vec))
                                    - 1/dot((Kink.j.vec-orb_b.vec),(Kink.j.vec-orb_b.vec)))
        delta_di -= delta_Tau * (lambda(e.N,e.rs)/2) *
                                (1/dot((Kink.k.vec-orb_a.vec),(Kink.k.vec-orb_a.vec))
                                    - 1/dot((Kink.k.vec-orb_b.vec),(Kink.k.vec-orb_b.vec)))
        delta_di -= delta_Tau * (lambda(e.N,e.rs)/2) *
                                (1/dot((Kink.l.vec-orb_a.vec),(Kink.l.vec-orb_a.vec))
                                    - 1/dot((Kink.l.vec-orb_b.vec),(Kink.l.vec-orb_b.vec)))
        Kink_semi_token = advance((c.kinks,Kink_semi_token))
        if Kink_semi_token == pastendsemitoken(c.kinks)
            Kink_semi_token = startof(c.kinks)
        end
        Tau_Kink,Kink = deref((c.kinks,Kink_semi_token))
        loop_counter += 1
    end
    return -delta_di * e.beta
end
