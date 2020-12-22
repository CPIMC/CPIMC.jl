function move_particle(c::Configuration, e::Ensemble)

    x = rand(c.occupations)
    oe = setdiff!(get_sphere_with_same_spin(x), c.occupations)
    if length(oe) == 0
        return 1
    end

    y = rand(oe)
    @assert x != y "same Configuration proposed."

    # weight change
    dw = exp(-e.β*(get_energy(y)-get_energy(x)))

    # change Configuration
    delete!(c.occupations,x)
    push!(c.occupations,y)

    # get orbitals for reverse update
    oe2 = setdiff!(get_sphere_with_same_spin(y), c.occupations)

    # quotient of proposal probabilities
    dv = length(oe)/length(oe2)

    dv*dw
end

function add_particle()
end

function remove_particle()
end
