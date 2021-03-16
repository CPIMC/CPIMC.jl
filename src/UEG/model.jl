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

OrbitalHEG(v::Tuple,s=0) = OrbitalHEG(SVector(v),s)

" returns a tuple of number of particles with spin up and number of particles with spin down "
function get_spin_up_down_count(c::Configuration)
    up = 0
    down = 0
    for orb in c.occupations
        if orb.spin == 1
            up += 1
        elseif orb.spin == -1
            down += 1
        else
            @assert(false)#Spin is not 1 or -1
        end
    end
    return (up,down)
end

function get_β_internal(θ::Float64, N::Int, c::Union{Nothing,Configuration} = nothing )
    # assuming that θ is defined with the Fermi energy corresponding to the particle number of the stronger occupied spin state
    # This is the case if the system is either fully spin-polarized or unpolarized
    if isnothing(c) # no configuration provided
        Spin_Faktor = 1
    else
        Spin_Faktor = max(get_spin_up_down_count(c)...)/N
    end
    return ((2*pi)^2)/(((6*(pi^2)*N*Spin_Faktor)^(2/3))*θ)
end

function λ(N::Int, rs::Float64)
    # TODO: why factor 2
    return (4/((2*pi)^3)) * (((4*pi)/3)^(1/3)) * rs * N^(1/3) * 2
end

function get_offdiagonal_element(e::Ensemble,kink::T4{OrbitalHEG{3}})
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
    # the factor λ/2 is due to the use of internal units
    wijkl *= λ(e.N,e.rs)/2
    # We sample with the weight of antisymmetrized matrix element but we do not restrict
    # the order of indizies of our possible kinks. We therefor need an extra factor 1/4 in the weight-function
    return wijkl * 1/4
end

function get_abs_offdiagonal_element(e::Ensemble,kink::T4{OrbitalHEG{3}})
    return abs(get_offdiagonal_element(e::Ensemble,kink::T4{OrbitalHEG{3}}))
end

" return the energy of a single diagonal single-particle matrix element "
function get_energy(o::OrbitalHEG)
    dot(o.vec,o.vec)
end

" coulomb kernel for 3D plane wavevectors "
kernel(i::StaticVector{N,Int}, k::StaticVector{N,Int}) where {N} = 1.0 / dot( i-k, i-k )

" coulomb kernel in plane wave basis "
kernel(i::OrbitalHEG,k::OrbitalHEG) = kernel(i.vec,k.vec)

## method declarations for the off-diagonal energy given by kinks in the UEG
# TODO: It may not be good style to use the same function name for the off-diagnoal energy as for the diagonal energy unless we dispatch at some point.
#       On the other hand, this notation is semantically convenient since it is short and conceptually unambiguous,
#       i.e. it is clear what is the difference between the energy of an orbital (diagonal) and the energy of a kink (off-diagonal).

""" return the energy of a single off-diagonal single-particle matrix element (a.k.a. kink)
    for a two-particle excitation (i.e. T4)"""
function get_energy(κ::T4)
    @assert iszero(κ.i + κ.j + κ.k + κ.l) "Momentum conservation violated by kink $(κ.k),$(κ.l) → $(κ.i), $(κ.j). All excitations must conserve the total momentum in the UEG."
    @assert !iszero(κ.i - κ.k) "Divergence in off-diagonal single-particle matrix element. All excitations must conserve the total momentum in the UEG."
    if κ.i.spin == κ.j.spin
        abs( kernel(κ.i, κ.k) - kernel(κ.i, κ.l) )# TODO absolute ?
    elseif κ.i.spin == κ.k.spin
        kernel(κ.i, κ.k)
    elseif κ.i.spin == κ.l.spin
        kernel(κ.i, κ.l)
    else
        @assert false "wrong spin combination in kink $(κ.k),$(κ.l) → $(κ.i), $(κ.j)."
    end
end

""" return the energy of a single off-diagonal matrix element (a.k.a. kink)
    for a one-particle excitation (i.e. T2)"""
function get_energy(κ::T2)
    throw(ErrorException("There are no T2-kinks in the UEG."))
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

" return two spheres of each spin with radius dk around the wavevector of the argument orbital o"
get_sphere(o::OrbitalHEG; dk=2) = union(get_sphere_with_same_spin(o,dk=dk),get_sphere_with_same_spin(OrbitalHEG(o.vec,Int(iszero(o.spin))),dk=dk))

function get_orbs_with_spin(orbitals::Set{OrbitalHEG{3}},spin::Int)
    orbs_s = Set{OrbitalHEG{3}}()
    for orb in orbitals
        if orb.spin==spin
            push!(orbs_s, orb)
        end
    end
    return orbs_s
end


""" calculates the change in the diagonal interaction
    when changing the ocupations between τ1 and τ2 accoring to left_kink
    already multiplied by e.β """
function get_change_diagonal_interaction(c::Configuration, e::Ensemble, left_kink::T4, τ1, τ2)
    orb_a = left_kink.i
    orb_b = left_kink.j
    orb_c = left_kink.k
    orb_d = left_kink.l
    delta_τ12 = τ2 - τ1
    if delta_τ12 < 0
        delta_τ12 += 1
    end
    if orb_a.spin == orb_b.spin
        delta_di = delta_τ12 * (λ(e.N,e.rs)/2) * (1/dot((orb_a.vec-orb_b.vec),(orb_a.vec-orb_b.vec)) -
                                            1/dot((orb_c.vec-orb_d.vec),(orb_c.vec-orb_d.vec)))
    else
        delta_di = 0
    end
    occs = occupations(c, τ1)
    for occ in occs
        if ((occ.vec == orb_c.vec) | (occ.vec == orb_d.vec) | (occ.vec == orb_a.vec) | (occ.vec == orb_b.vec))
            nothing
        else
            for orb in [orb_a, orb_b]
                if occ.spin == orb.spin
                    delta_di += delta_τ12 * (λ(e.N,e.rs)/2) * 1/dot((occ.vec-orb.vec),(occ.vec-orb.vec))
                end
            end
            for orb in [orb_c, orb_d]
                if occ.spin == orb.spin
                    delta_di -= delta_τ12 * (λ(e.N,e.rs)/2) * 1/dot((occ.vec-orb.vec),(occ.vec-orb.vec))
                end
            end
        end
        @assert(abs(delta_di) != Inf)
    end
    if length(c.kinks) == 0
        return -delta_di * e.β
    end
    kink_semi_token = searchsortedfirst(c.kinks,τ1)
    if kink_semi_token == pastendsemitoken(c.kinks)
        kink_semi_token = startof(c.kinks)
    end
    τ_kink,kink = deref((c.kinks,kink_semi_token))
    # the kink at τ1 is already considered in occs
    if τ_kink == τ1
        kink_semi_token = advance((c.kinks,kink_semi_token))
        if kink_semi_token == pastendsemitoken(c.kinks)
            kink_semi_token = startof(c.kinks)
        end
        τ_kink,kink = deref((c.kinks,kink_semi_token))
    end
    loop_counter = 0
    @assert(abs(delta_di) != Inf)
    while ((τ1 < τ_kink < τ2) | (τ_kink < τ2 < τ1) | (τ2 < τ1 < τ_kink)) & (loop_counter < length(c.kinks))
        delta_τ = τ2 - τ_kink
        if delta_τ < 0
            delta_τ += 1
        end
        for occ in [kink.i, kink.j]
            for orb in [orb_a, orb_b]
                if occ.spin == orb.spin
                    delta_di += delta_τ * (λ(e.N,e.rs)/2) *
                                            (1/dot((occ.vec-orb.vec),(occ.vec-orb.vec)))
                end
            end
        end
        for occ in [kink.k, kink.l]
            for orb in [orb_c, orb_d]
                if occ.spin == orb.spin
                    delta_di += delta_τ * (λ(e.N,e.rs)/2) *
                                            (1/dot((occ.vec-orb.vec),(occ.vec-orb.vec)))
                end
            end
        end
        for occ in [kink.i, kink.j]
            for orb in [orb_c, orb_d]
                if occ.spin == orb.spin
                    delta_di -= delta_τ * (λ(e.N,e.rs)/2) *
                                            (1/dot((occ.vec-orb.vec),(occ.vec-orb.vec)))
                end
            end
        end
        for occ in [kink.k, kink.l]
            for orb in [orb_a, orb_b]
                if occ.spin == orb.spin
                    delta_di -= delta_τ * (λ(e.N,e.rs)/2) *
                                            (1/dot((occ.vec-orb.vec),(occ.vec-orb.vec)))
                end
            end
        end
        @assert(abs(delta_di) != Inf)
        kink_semi_token = advance((c.kinks,kink_semi_token))
        if kink_semi_token == pastendsemitoken(c.kinks)
            kink_semi_token = startof(c.kinks)
        end
        τ_kink,kink = deref((c.kinks,kink_semi_token))
        loop_counter += 1
    end
    return -delta_di * e.β
end

" This function assumes for now that all kinks in configuration.kinks are of type 4. "
function get_change_diagonal_interaction(c::Configuration, e::Ensemble, left_kink::T2, τ1, τ2)
    orb_a = left_kink.i
    orb_b = left_kink.j
    delta_τ12 = τ2 - τ1
    if delta_τ12 < 0
        delta_τ12 += 1
    end
    delta_di = 0
    occs = occupations(c, τ1)
    @assert(!in(orb_a, occs))
    for occ in occs
        if ((occ.vec == orb_b.vec) | (occ.vec == orb_a.vec))
            nothing
        else
            if occ.spin == orb_a.spin
                delta_di += delta_τ12 * (λ(e.N,e.rs)/2) * 1/dot((occ.vec-orb_a.vec),(occ.vec-orb_a.vec))
            end
            if occ.spin == orb_b.spin
                delta_di -= delta_τ12 * (λ(e.N,e.rs)/2) * 1/dot((occ.vec-orb_b.vec),(occ.vec-orb_b.vec))
            end
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
    # The kink at τ1 is already considered in occs.
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

        for occ in [kink.i, kink.j]
            if occ.spin == orb_a.spin
                delta_di += delta_τ * (λ(e.N,e.rs)/2) *
                                        (1/dot((occ.vec-orb_a.vec),(occ.vec-orb_a.vec)))
            end
        end
        for occ in [kink.k, kink.l]
            if occ.spin == orb_b.spin
                delta_di += delta_τ * (λ(e.N,e.rs)/2) *
                                        (1/dot((occ.vec-orb_b.vec),(occ.vec-orb_b.vec)))
            end
        end
        for occ in [kink.i, kink.j]
            if occ.spin == orb_b.spin
                delta_di -= delta_τ * (λ(e.N,e.rs)/2) *
                                        (1/dot((occ.vec-orb_b.vec),(occ.vec-orb_b.vec)))
            end
        end
        for occ in [kink.k, kink.l]
            if occ.spin == orb_a.spin
                delta_di -= delta_τ * (λ(e.N,e.rs)/2) *
                                        (1/dot((occ.vec-orb_a.vec),(occ.vec-orb_a.vec)))
            end
        end

        kink_semi_token = advance((c.kinks,kink_semi_token))
        if kink_semi_token == pastendsemitoken(c.kinks)
            kink_semi_token = startof(c.kinks)
        end
        τ_kink,kink = deref((c.kinks,kink_semi_token))
        loop_counter += 1
    end
    return -delta_di * e.β
end

function get_sign_offdiagonal_product(e,c)
    sign_ofd_prod = 1
    for (_,kink) in c.kinks
        sign_ofd_prod *= sign(get_offdiagonal_element(e,kink))
    end
    return sign_ofd_prod
end
