using CPIMC
using CPIMC.PlaneWaves

import CPIMC: energy, Step
import CPIMC.UniformElectronGas: λ, β

import LinearAlgebra: dot

struct FreeElectronGas <: Model end

energy(::FreeElectronGas, o::PlaneWave) = dot(o.vec, o.vec)

function move_particle(m::Model, e::Ensemble, c::Configuration)
    x = rand(c.occupations)
    oe = setdiff!(sphere_with_same_spin(x), c.occupations)

    if isempty(oe)
        return 1.0, Step()
    end
    y = rand(oe)

    # weight change
    dw = exp(-e.β*(energy(m,y)-energy(m,x)))

    # get orbitals for reverse update
    oe2 = setdiff!(sphere_with_same_spin(y), c.occupations)

    # quotient of proposal probabilities
    δv = length(oe)/length(oe2)

    δv * dw, Step(x,y)
end

function Ekin(m::Model, e::Ensemble, c::Configuration) :: UInt
    sum(energy(m, n) for n in c.occupations)
end

function occupations(m::Model, e::Ensemble, c::Configuration, emax::Int=100) :: Array{UInt,1}
    nk = zeros(UInt, emax)

    ens = [ energy(m, n) for n in c.occupations ]

    for ε in ens[ens .< emax]
        nk[ε+1] = nk[ε+1] + 1
    end
    nk
end


function main()
    # MC options
    NMC = 10^6
    cyc = 10
    NEquil = 10^4
    # system parameters
    θ = 1.0
    rs = 0.5

    S = sphere_with_same_spin(PlaneWave((0,0,0),Up),dk=2)### use 33 particles
    N = length(S)
    ξ = fractional_spin_polarization(S)
    c = Configuration(S)
    d = dimension(c.occupations)

    println("θ: ", θ)
    println("rs: ", rs)
    println("N: ", N)
    println("d: ", d)

    e = CEnsemble(λ(N,rs,d), β(θ,N,ξ,d), N) # β only works for 3D

    update_names = [move_particle]
    updates = Array{Tuple{Function,UpdateCounter},1}()
    for up_name in update_names
        push!(updates, (up_name,UpdateCounter()))
    end

    measurements = Dict(
      :Ekin => (Variance(UInt), Ekin)
    , :occs => (Group([Variance(UInt) for i in 1:100]), occupations))

    println("Start MC process ... ")
    sweep!(FreeElectronGas(), e, c, updates, measurements, NMC, cyc, NEquil)
    println(" finished.")

    print_results(measurements, e)

end

main()
