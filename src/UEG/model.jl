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



# the remark on the field may be removed if it considered to be clear that OrbitalHEG must always have a field spin, or if a function `spin(o::Orbital)` is used
"""
    fractional_spin_polarization(occ::Set{OrbitalHEG)
calculate the fractional spin polarization from a set of Orbitals which each must have a field `spin` of type `@enum Spin Down Up`,
the fractional spin polarization is defined as the ratio of the absolute difference of the number of particles which Spin Down and Spin Up and the particle number,
thus here no convention is made as to which spin component should be occupied more often, as opposed to some literature where N↑ > N↓ is used
"""
function fractional_spin_polarization(occ::Set{<:OrbitalHEG})
    N_up = length(filter(x -> x.spin == Up, occ))
    N_down = length(filter(x -> x.spin == Down, occ))
    N = length(occ)
    return abs(N_up - N_down) / N
end

"""
    ξ(occ::Set{OrbitalHEG})
calculate the fractional spin polarization from a set of Orbitals which each must have a field `spin` of type `@enum Spin Down Up`,
the fractional spin polarization is defined as the ratio of the absolute difference of the number of particles which Spin Down and Spin Up and the particle number,
thus here no convention is made as to which spin component should be occupied more often, as opposed to some literature where N↑ > N↓ is used
"""
ξ(occ::Set{<:OrbitalHEG}) = fractional_spin_polarization(occ)


"""
    fractional_spin_polarization(c::Configuration{<:OrbitalHEG})
calculate the fractional spin polarization from a set of Orbitals which each must have a field `spin` of type `@enum Spin Down Up`,
the fractional spin polarization is defined as the ratio of the absolute difference of the number of particles which Spin Down and Spin Up and the particle number,
thus here no convention is made as to which spin component should be occupied more often, as opposed to some literature where N↑ > N↓ is used
"""
function fractional_spin_polarization(c::Configuration{<:OrbitalHEG})# TODO: added {<:OrbitalHEG}
    fractional_spin_polarization(c.occupations)
end


"""
    ξ(c::Configuration{<:OrbitalHEG})
calculate the fractional spin polarization from a set of Orbitals which each must have a field `spin` of type `@enum Spin Down Up`,
the fractional spin polarization is defined as the ratio of the absolute difference of the number of particles which Spin Down and Spin Up and the particle number,
thus here no convention is made as to which spin component should be occupied more often, as opposed to some literature where N↑ > N↓ is used
"""
ξ(c::Configuration{<:OrbitalHEG}) = fractional_spin_polarization(c)


"""
    λ(N, rs, d=3)

calculate λ, which is the coupling constant of the Hamiltonian in internal units
from the particle number N, Bruecker parameter rs for the density and spatial dimension d
the differences for d ∈ {2,3} reflect the different expressions for the Coulomb interaction, whichs Fourier component is v_q = \frac{4π}{q^2} in 3D and v_q = \frac{2π}{q} in 2D
the one-dimensional case involves a short-distance cutoff to avoid the non-integrable divergence of the interaction at zero distance and is not implemented
"""
# The dimension is given as an optional argument with default value 3 in order to work with the current implementation, which requires λ to be defined for 2 arguments (N, rs). When λ is part of struct Ensemble, the default value can be removed.
function λ(N, rs, d::Int=3)
    if d == 3
        (4 / (2π)^3 ) * cbrt(4π*N/3) * rs
    elseif d == 2
        sqrt(N/π) * rs / (2π)
    else
        throw(DomainError("coupling constant λ is only implemented for dimension d ∈ {2,3} and not for d=$(d)"))
    end
end


"""
    rs(N::Int, λ::Float64, d=3)
calculate the Brueckner parameter rs for the density from the particle number N
and λ, which is the coupling constant of the Hamiltonian in internal units"""
# The dimension is given as an optional argument with default value 3 in order to work with the current implementation, which requires λ to be defined for 2 arguments (N, rs). When λ is part of struct Ensemble, the default value can be removed.
function rs(N, λ, d::Int=3)
    if d == 3
        ( (2π)^3 / 4 ) * cbrt(3 / (4π*N)) * λ
    elseif d == 2
        2π * sqrt(π/N) * λ
    else
        throw(DomainError("Coupling constant λ is only implemented for dimension d ∈ {2,3} and not for d=$(d). Thus, no conversion from λ to rs is defined for d ∉ {2,3}"))
    end
end


"""
    internal_energy_factor(N, rs, d=3)

This factor is used to convert an energy quantity from Hartree a.u. to internal units.
A value of the quantity in Hartree a.u. has to be multiplied by this factor to yield the corresponding value of the quantity in internal units.
A value of the quantity in internal units has to be divided by this factor to yield the corresponding value of the quantity in Hartree a.u.
"""
# The dimension is given as an optional argument with default value 3 in order to work with the current implementation, which requires λ to be defined for 2 arguments (N, rs). When λ is part of struct Ensemble, the default value can be removed.
internal_energy_factor(N, rs, d=3) = (2π)^4 * λ(N, rs, d)^2 / 8

"""
    kF(rs, ξ, d::Int=3)
calculate the Fermi wavenumber in Hartree a.u. from Brueckner parameter rs for the density,
the fractional spin polarization ξ and spatial dimension d
"""
# The dimension is given as an optional argument with default value 3 in order to work with the current implementation, which requires λ to be defined for 2 arguments (N, rs). When λ is part of struct Ensemble, the default value can be removed.
function kF(rs, ξ, d::Int=3)
    if d == 1
        α = 4 / π
    elseif d == 2
        α = 1 / sqrt(2)
    elseif d == 3
        α = cbrt( 4 / (9π) )
    else
        throw(DomainError("Fermi-wavenumber only defined for 1-, 2- and 3-dimensional systems. Choose d ∈ {1,2,3} and not d=$(d)"))
    end
    (1 + ξ)^(1/d) / ( α*rs )
end


"""
    EF(rs, ξ, d::Int=3)
calculate the Fermi energy in Hartree a.u. from Brueckner parameter rs for the density,
the fractional spin polarization ξ and spatial dimension d
"""
# The dimension is given as an optional argument with default value 3 in order to work with the current implementation, which requires λ to be defined for 2 arguments (N, rs). When λ is part of struct Ensemble, the default value can be removed.
EF(rs, ξ, d::Int=3) = kF(rs, ξ, d)^2 / 2

"""
    β_au(Θ, EF)
calculate the inverse temperature in Hartree a.u. from the reduced temperature Θ and the Fermi energy EF in Hartree a.u.
the inverse temperature is defined by β = 1 / (kB*T) where kB is the Boltzmann constant and T is the temperature
"""
β_Ha(Θ, EF) = 1 / (Θ * EF)

"""
    β(Θ, rs, N, ξ, d=3)

calculate β in internal units
"""
# The dimension is given as an optional argument with default value 3 in order to work with the current implementation, which requires λ to be defined for 2 arguments (N, rs). When λ is part of struct Ensemble, the default value can be removed.
β(Θ, rs, N, ξ, d=3) = 1 / ( Θ * EF(rs, ξ, d) * internal_energy_factor(N, rs, d) )


" coulomb kernel for 3D plane wavevectors "
kernel(i::StaticVector{N,Int}, k::StaticVector{N,Int}) where {N} = 1.0 / dot( i-k, i-k )

" coulomb kernel in plane wave basis "
kernel(i::OrbitalHEG,k::OrbitalHEG) = kernel(i.vec,k.vec)




" anti-symmetrized interaction matrix element "
function wminus(i::OrbitalHEG{D}, j::OrbitalHEG{D}, k::OrbitalHEG{D}, l::OrbitalHEG{D}) where {D}
    @assert ((i != k) && (i != l))
    @assert ((i.spin == j.spin) == (k.spin == l.spin))
    @assert(in(i.spin,[k.spin,l.spin]))
    @assert(i.vec + j.vec == k.vec + l.vec)
    if i.spin == j.spin
        @assert(!isinf(abs(kernel(i, k) - kernel(i, l))))
        return kernel(i, k) - kernel(i, l)
    elseif i.spin == k.spin
        @assert(!isinf(abs(kernel(i, k))))
        return kernel(i, k)
    elseif i.spin == l.spin
        @assert(!isinf(abs(kernel(i, l))))
        return -kernel(i, l)
    end
end

" diagonal interaction matrix element "
function wdiag(a::OrbitalHEG{D}, b::OrbitalHEG{D}) where {D}
    if a.spin != b.spin
        return 0
    else
        return -kernel(a, b)
    end
end

# " anti-symmetrized interaction matrix element "
# function wminus(i::OrbitalHEG{D}, j::OrbitalHEG{D}, k::OrbitalHEG{D}, l::OrbitalHEG{D}) where {D}
#     if (i.vec == k.vec)
#         @assert j.vec == l.vec
#         return 0.0
#     elseif (i.vec == l.vec)
#         @assert(j.vec == k.vec)
#         if i.spin == j.spin
#             @assert(k.spin == l.spin)
#             return -kernel(i,k)
#         else
#             return 0.0
#         end
#     else
#         if i.spin == j.spin
#             return kernel(i, k) - kernel(i, l)
#         elseif i.spin == k.spin
#             return kernel(i, k)
#         elseif i.spin == l.spin
#             return -kernel(i, l)
#         end
#     end
# end

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

    Δdi = Δτ12 * λ(e.N,e.rs) * ( wdiag(orb_a, orb_b) - wdiag(orb_c, orb_d) )

    occs = occupations(c, τ1)

    # collect diagonal interaction energy at τ1
    for occ in occs
        if occ.vec in [ orb_a.vec, orb_b.vec, orb_c.vec, orb_d.vec ]
            continue
        else
            for orb in [orb_a, orb_b]
                Δdi += Δτ12 * λ(e.N,e.rs) * wdiag(occ,orb)
            end
            for orb in [orb_c, orb_d]
                Δdi -= Δτ12 * λ(e.N,e.rs) * wdiag(occ,orb)
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
                 Δdi += Δτ * λ(e.N,e.rs) * wdiag(occ,orb)
            end
            for orb in [orb_c, orb_d]
                 Δdi -= Δτ * λ(e.N,e.rs) * wdiag(occ,orb)
            end
        end
        for occ in [kink.k, kink.l]
            for orb in [orb_a, orb_b]
                 Δdi -= Δτ * λ(e.N,e.rs) * wdiag(occ,orb)
            end
            for orb in [orb_c, orb_d]
                 Δdi += Δτ * λ(e.N,e.rs) * wdiag(occ,orb)
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
            Δdi += Δτ12 * λ(e.N,e.rs) * wdiag(occ,orb_a)
            Δdi -= Δτ12 * λ(e.N,e.rs) * wdiag(occ,orb_b)
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
            Δdi += Δτ * λ(e.N,e.rs) * wdiag(occ,orb_a)
            Δdi -= Δτ * λ(e.N,e.rs) * wdiag(occ,orb_b)
        end
        for occ in [kink.k, kink.l]
            Δdi += Δτ * λ(e.N,e.rs) * wdiag(occ,orb_b)
            Δdi -= Δτ * λ(e.N,e.rs) * wdiag(occ,orb_a)
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
