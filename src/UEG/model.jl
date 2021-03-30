using StaticArrays

import LinearAlgebra: dot

include("orbital.jl")

#TODO: store λ instead of rs
struct Ensemble_UEG <: Ensemble
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
function β(Θ, N, ξ, d=3)
    #calculate the factor alpha used in calculation of the fermi-vector in Ha units
    if d == 1
        α = 4 / π
    elseif d == 2
        α = 1 / sqrt(2)
    elseif d == 3
        α = cbrt( 4 / (9π) )
    else
        throw(DomainError("Fermi-wavenumber only defined for 1-, 2- and 3-dimensional systems. Choose d ∈ {1,2,3} and not d=$(d)"))
    end

    #calculate λ_rs_ratio
    if d == 3
        λ_rs_ratio = (4 / (2π)^3 ) * cbrt(4π*N/3)
    elseif d == 2
        λ_rs_ratio = sqrt(N/π) / (2π)
    else
        throw(DomainError("coupling constant λ is only implemented for dimension d ∈ {2,3} and not for d=$(d)"))
    end

    return (16/Θ) * (α / ((1 + ξ)^(1/d) * (2π)^2 * λ_rs_ratio))^2
end


" coulomb kernel for 3D plane wavevectors "
kernel(i::StaticVector{N,Int}, k::StaticVector{N,Int}) where {N} = 1.0 / dot( i-k, i-k )

" coulomb kernel in plane wave basis "
kernel(i::OrbitalHEG,k::OrbitalHEG) = kernel(i.vec,k.vec)
