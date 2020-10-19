function Ekin(c::Configuration) :: Float64
    sum(get_energy(n) for n in c.occupations)
end

function occupations(c::Configuration, emax::Int=100) :: Array{UInt,1}
    nk = zeros(UInt, emax)

    ens = get_energy.(c.occupations)

    for ε in ens[ens .< emax]
        nk[ε+1] = nk[ε+1] + 1
    end
    nk
end

function particleNumber(c::Configuration) :: Int
  return length(c.occupations)
end

abstract type Observable end
abstract type Energy <: Observable end
abstract type Occupation <: Observable end

abstract type UnitSystem end
abstract type InternalUnits <: UnitSystem end
abstract type HartreeUnits <: UnitSystem end
abstract type SIUnits <: UnitSystem end

abstract type Quantity{O<:Observable,U<:UnitSystem} end
