using StaticArrays

import LinearAlgebra: dot

include("orbital.jl")

#TODO: store λ instead of rs
struct Ensemble
  "Brueckner parameter"
  rs :: Float64
  "reduced temperature"
  β :: Float64
  "particle number"
  N :: Int
end


function fractional_spin_polarization(occ::Set{OrbitalHEG{D}}) where {D}
    up = 0
    down = 0

    for orb in occ
        if orb.spin == Up
            up += 1
        else
            down += 1
        end
    end

    return (up-down)/(up+down)
end

function fractional_spin_polarization(c::Configuration)
    fractional_spin_polarization(c.occupations)
end


"""
    β(θ::Float64, N::Int, ξ::Float64 )

calculate β in internal units
"""
function β(θ::Float64, N::Int, ξ::Float64 = 1.0)
    return (2*pi)^2/(((6*(pi^2) * N/2 * (1+abs(ξ)))^(2/3))*θ)
end

"""
    λ(N::Int, rs::Float64)

calculate λ
"""
function λ(N::Int, rs::Float64)
    return 4/((2*pi)^3) * (4*pi/3)^(1/3) * rs * N^(1/3)
end


" coulomb kernel for 3D plane wavevectors "
kernel(i::StaticVector{N,Int}, k::StaticVector{N,Int}) where {N} = 1.0 / dot( i-k, i-k )

" coulomb kernel in plane wave basis "
kernel(i::OrbitalHEG,k::OrbitalHEG) = kernel(i.vec,k.vec)


" anti-symmetrized interaction matrix element "
function wminus(i::OrbitalHEG{D}, j::OrbitalHEG{D}, k::OrbitalHEG{D}, l::OrbitalHEG{D}) where {D}
    if (i.vec == k.vec)
        @assert j.vec == l.vec)
        return 0.0
    elseif (i.vec == l.vec)
        @assert(j.vec == k.vec)
        if i.spin == j.spin
            @assert(k.spin == l.spin)
            return -kernel(i,k)
        else
            return 0.0
        end
    else
        if i.spin == j.spin
            return kernel(i, k) - kernel(i, l)
        elseif i.spin == k.spin
            return kernel(i, k)
        elseif i.spin == l.spin
            return -kernel(i, l)
        end
    end
end

wminus(kink::T4) = wminus(kink.i, kink.j, kink.k, kink.l)

function offdiagonal_element(e::Ensemble, kink::T4)
    # We sample with the weight of antisymmetrized matrix element but we do not restrict
    # the order of indices of our possible kinks. We therefor need an extra factor 1/4 in the weight-function
    return 1/4 * λ(e.N,e.rs) * wminus(kink)
end



"""
    Δdiagonal_interaction(c::Configuration, e::Ensemble, left_kink::T4, τ1, τ2)

change in the diagonal part of the interaction when applying left_kink between τ1 and τ2
"""

function Δdiagonal_interaction(c::Configuration, e::Ensemble, orb_a::Orbital, orb_b::Orbital, orb_c::Orbital, orb_d::Orbital, τ1, τ2)
    Δτ12 = τ2 - τ1

    if Δτ12 < 0
        Δτ12 += 1
    end

    @assert (orb_a.spin != orb_b.spin) == (orb_c.spin != orb_d.spin )

    Δdi = Δτ12 * λ(e.N,e.rs) * ( wminus(orb_a, orb_b, orb_b, orb_a) - wminus(orb_c, orb_d, orb_d, orb_c) )

    occs = occupations(c, τ1)

    # collect diagonal interaction energy at τ1
    for occ in occs
        if occ.vec in [ orb_a.vec, orb_b.vec, orb_c.vec, orb_d.vec ]
            continue
        else
            for orb in [orb_a, orb_b]
                Δdi += Δτ12 * λ(e.N,e.rs) * wminus(occ,orb,orb,occ)
            end
            for orb in [orb_c, orb_d]
                Δdi -= Δτ12 * λ(e.N,e.rs) * wminus(occ,orb,orb,occ)
            end
        end
        @assert !isinf(abs(Δdi))
        @assert(!isnan(Δdi))
    end

    if isempty(c.kinks)
        return Δdi
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

    @assert !isinf(abs(Δdi))
    @assert(!isnan(Δdi))
    # collect contributions to diagonal interaction energy due to kinks in the interval
    while ((τ1 < τ_kink < τ2) | (τ_kink < τ2 < τ1) | (τ2 < τ1 < τ_kink)) & (loop_counter < length(c.kinks))
        Δτ = τ2 - τ_kink
        if Δτ < 0
            Δτ += 1
        end
        for occ in [kink.i, kink.j]
            for orb in [orb_a, orb_b]
                 Δdi += Δτ * λ(e.N,e.rs) * wminus(occ,orb,orb,occ)
            end
            for orb in [orb_c, orb_d]
                 Δdi -= Δτ * λ(e.N,e.rs) * wminus(occ,orb,orb,occ)
            end
        end
        for occ in [kink.k, kink.l]
            for orb in [orb_a, orb_b]
                 Δdi -= Δτ * λ(e.N,e.rs) * wminus(occ,orb,orb,occ)
            end
            for orb in [orb_c, orb_d]
                 Δdi += Δτ * λ(e.N,e.rs) * wminus(occ,orb,orb,occ)
            end
        end

        @assert !isinf(abs(Δdi))
        @assert(!isnan(Δdi))
        kink_semi_token = advance((c.kinks,kink_semi_token))
        if kink_semi_token == pastendsemitoken(c.kinks)
            kink_semi_token = startof(c.kinks)
        end
        τ_kink,kink = deref((c.kinks,kink_semi_token))
        loop_counter += 1
    end
    @assert(!isnan(Δdi))
    return Δdi
end

" This function assumes for now that all kinks in c are of type 4. "
function Δdiagonal_interaction(c::Configuration, e::Ensemble, orb_a::Orbital, orb_b::Orbital, τ1, τ2)
    Δτ12 = τ2 - τ1
    if Δτ12 < 0
        Δτ12 += 1
    end
    Δdi = 0
    occs = occupations(c, τ1)
    @assert !in(orb_a, occs)
    for occ in occs
        if occ.vec in [ orb_a.vec, orb_b.vec ]
            continue
        else
            Δdi += Δτ12 * λ(e.N,e.rs) * wminus(occ,orb_a,orb_a,occ)
            Δdi -= Δτ12 * λ(e.N,e.rs) * wminus(occ,orb_b,orb_b,occ)
        end
    end
    if isempty(c.kinks)
        return Δdi
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
        Δτ = τ2 - τ_kink
        if Δτ < 0
            Δτ += 1
        end

        for occ in [kink.i, kink.j]
            Δdi += Δτ * λ(e.N,e.rs) * wminus(occ,orb_a,orb_a,occ)
            Δdi -= Δτ * λ(e.N,e.rs) * wminus(occ,orb_b,orb_b,occ)
        end
        for occ in [kink.k, kink.l]
            Δdi += Δτ * λ(e.N,e.rs) * wminus(occ,orb_b,orb_b,occ)
            Δdi -= Δτ * λ(e.N,e.rs) * wminus(occ,orb_a,orb_a,occ)
        end


        kink_semi_token = advance((c.kinks,kink_semi_token))
        if kink_semi_token == pastendsemitoken(c.kinks)
            kink_semi_token = startof(c.kinks)
        end
        τ_kink,kink = deref((c.kinks,kink_semi_token))
        loop_counter += 1
    end

    return Δdi
end

"""
    sign_offdiagonal_product(c::Configuration)

sign of the product of the offdiagonal elements
used for the calculation the sign of the weight function
"""
function sign_offdiagonal_product(c::Configuration)
    sign_ofd_prod = 1
    for kink in values(c.kinks)
        sign_ofd_prod *= sign(wminus(kink))
    end
    return sign_ofd_prod
end
