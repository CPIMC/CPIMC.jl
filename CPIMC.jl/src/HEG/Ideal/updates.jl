function move_particle(c::Configuration, e::Ensemble)

    x = rand(c.occupations)
    oe = get_orbs_with_spin(setdiff!(get_sphere(x), c.occupations),x.spin)
    #@assert length(oe) > 0 "no empty orbitals in neighbourhood"
    if length(oe) == 0
        return 1
    end

    y = rand(oe)
    @assert x != y "same Configuration proposed."

    # weight change
    dw = exp(-e.beta*(get_energy(y)-get_energy(x)))

    # change Configuration
    delete!(c.occupations,x)
    push!(c.occupations,y)

    # get orbitals for reverse update
    oe2 = get_orbs_with_spin(setdiff!(get_sphere(y), c.occupations),y.spin)

    # quotient of proposal probabilities
    dv = length(oe)/length(oe2)

    dv*dw
end

function add_particle()
end

function remove_particle()
end
