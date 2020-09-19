function move_particle(c::Configuration, e::Ensemble)

    x = rand(c.occupations)
    oe = setdiff!(get_sphere(x), c.occupations)
    @assert length(oe) > 0 "no empty orbitals in neighbourhood"

    y = rand(oe)
    @assert x != y "same Configuration proposed."

    # weight change
    dw = exp(-e.beta*(get_energy(y)-get_energy(x)))

    # change Configuration
    delete!(c.occupations,x)
    push!(c.occupations,y)

    # get orbitals for reverse update
    oe2 = setdiff!(get_sphere(y), c.occupations)

    # quotient of proposal probabilities
    dv = length(oe)/length(oe2)

    dv*dw
end

function add_particle()
end

function remove_particle()
end

