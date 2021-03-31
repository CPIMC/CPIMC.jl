module ElectronGas

using CPIMC
using CPIMC.PlaneWaves

import CPIMC: kernel, energy


import StaticArrays: StaticVector
import LinearAlgebra: dot

export UEG, β, EF, λ, β_Ha, rs, kF


struct UEG <: Model end

" return the energy of a single diagonal single-particle matrix element "
function energy(m::UEG, o::PlaneWave)
    dot(o.vec,o.vec)
end

" coulomb kernel in plane wave basis "
kernel(m::UEG, i::PlaneWave, k::PlaneWave) = 1.0 / dot(i.vec - k.vec, i.vec - k.vec)

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
