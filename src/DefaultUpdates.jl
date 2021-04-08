"""
This module implements a canonical set of updates for models with momentum-conservation and only two-particle excitations (no `T2` kinks).
Currently limited to choosing `PlaneWave`s as basis functions.
"""
module DefaultUpdates

export ex_radius, isuseful

"""
    const ex_radius = 3

radius of the sphere of orbitals which are considered for excitations
"""
const ex_radius = 3 # TODO: find better solution
# ex_radius could be kwarg to each update function
# and passed to `sweep!` via anonymous function with local variable ex_radius

using ..CPIMC
using ..CPIMC.PlaneWaves

import ..CPIMC: τ_borders, Woffdiag_element, isunaffected_in_interval, ΔWdiag_element, isunaffected, occupations_at, drop, add, energy, random_shuffle, basis, τ_next_affecting, τ_prev_affecting

using LinearAlgebra

include("updates/type_a.jl")
include("updates/type_b.jl")
include("updates/type_c.jl")
include("updates/type_d.jl")
include("updates/type_e.jl")

include("updates/auxiliary.jl")

end

