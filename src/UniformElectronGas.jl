"""
Types and definitions used for simulations of the uniform electron gas.
"""
module UniformElectronGas

using ..CPIMC
using ..CPIMC.PlaneWaves

import ..CPIMC: energy, w


import StaticArrays: StaticVector
import LinearAlgebra: dot

export UEG, β, EF, λ, β_Ha, rs, kF

"""
Singleton for dispatching.

For other models this would be the place to store additional
parameters that might come up in the matrix elements.
"""
struct UEG <: Model end

@doc raw"""
    energy(m::UEG, o::PlaneWave)

The single-particle energy of a plane wave with momentum $\vec{k}$ is given by (Hartree atomic units):

$ \epsilon_k = \vec{k}^2 / 2. $

In the internal unit system we use integer $k$-values, which results in:

$ \epsilon_k = \vec{k}^2. $

"""
function energy(m::UEG, o::PlaneWave)
    dot(o.vec,o.vec)
end

"""
    w(i::Orbital, j::Orbital, k::Orbital, l::Orbital)

Return the two-particle matrix element of the Coulomb interaction.
Zero is returned if momentum or spin is not conserved, in consequence of the Bloch-theorem
and the spin kronecker-delta in the plane-spin-wave basis.

The singular contribution `(i.vec == k.vec)` is also removed, i.e. set to zero.
This property is not a result of the Bloch theorem for the two-particle Coulomb matrix
element in the plane wave basis but is added here to allow for straightforward summation over
all combination of orbitals without explictly excluding this component.

It is not part of the UEG many-body Hamiltonian due to the physical assumption that contributions from the
uniform background cancel with this diverging term in the thermodynamic limit.
"""
function w(m::UEG, i::Orbital, j::Orbital, k::Orbital, l::Orbital)
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
    λ(N, rs, d::Int)

calculate λ, which is the coupling constant of the Hamiltonian in internal units
from the particle number N, Bruecker parameter rs for the density and spatial dimension d
the differences for d ∈ {2,3} reflect the different expressions for the Coulomb interaction, whichs Fourier component is v_q = \frac{4π}{q^2} in 3D and v_q = \frac{2π}{q} in 2D
the one-dimensional case involves a short-distance cutoff to avoid the non-integrable divergence of the interaction at zero distance and is not implemented
"""
function λ(N, rs, d::Int)
    if d == 3
        (4 / (2π)^3 ) * cbrt(4π*N/3) * rs
    elseif d == 2
        sqrt(N/π) * rs / (2π)
    else
        throw(DomainError("coupling constant λ is only implemented for dimension d ∈ {2,3} and not for d=$(d)"))
    end
end


"""
    rs(N, λ, d::Int)
calculate the Brueckner parameter rs for the density from the particle number N
and λ, which is the coupling constant of the Hamiltonian in internal units"""
function rs(N, λ, d::Int)
    if d == 3
        ( (2π)^3 / 4 ) * cbrt(3 / (4π*N)) * λ
    elseif d == 2
        2π * sqrt(π/N) * λ
    else
        throw(DomainError("Coupling constant λ is only implemented for dimension d ∈ {2,3} and not for d=$(d). Thus, no conversion from λ to rs is defined for d ∉ {2,3}"))
    end
end


"""
    internal_energy_factor(N, rs, d)

This factor is used to convert an energy quantity from Hartree a.u. to internal units.
A value of the quantity in Hartree a.u. has to be multiplied by this factor to yield the corresponding value of the quantity in internal units.
A value of the quantity in internal units has to be divided by this factor to yield the corresponding value of the quantity in Hartree a.u.
"""
internal_energy_factor(N, rs, d) = (2π)^4 * λ(N, rs, d)^2 / 8

"""
    kF(rs, ξ, d::Int)
calculate the Fermi wavenumber in Hartree a.u. from Brueckner parameter rs for the density,
the fractional spin polarization ξ and spatial dimension d
"""
function kF(rs, ξ, d::Int)
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
    EF(rs, ξ, d::Int)
calculate the Fermi energy in Hartree a.u. from Brueckner parameter rs for the density,
the fractional spin polarization ξ and spatial dimension d
"""
EF(rs, ξ, d::Int) = kF(rs, ξ, d)^2 / 2

"""
    β_Ha(Θ, EF)
calculate the inverse temperature in Hartree a.u. from the reduced temperature Θ and the Fermi energy EF in Hartree a.u.
the inverse temperature is defined by β = 1 / (kB*T) where kB is the Boltzmann constant and T is the temperature
"""
β_Ha(Θ, EF) = 1 / (Θ * EF)

"""
    β(Θ, rs, N, ξ, d::Int)

calculate β in internal units
"""
function β(Θ, N, ξ, d::Int)
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


end
