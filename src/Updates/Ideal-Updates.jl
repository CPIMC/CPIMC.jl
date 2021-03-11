function move_particle(c::Configuration, e::Ensemble) :: Tuple{Float64, Step}
    x = rand(c.occupations)
    oe = setdiff!(get_sphere_with_same_spin(x), c.occupations)
    if isempty(oe)
        return 1.0, Step()
    end
    y = rand(oe)
    @assert x != y "same Configuration proposed."

    # weight change
    dw = exp(-e.β*(get_energy(y)-get_energy(x)))

    # get orbitals for reverse update
    oe2 = setdiff!(get_sphere_with_same_spin(y), c.occupations)

    # quotient of proposal probabilities
    δv = length(oe)/length(oe2)

    δv * dw, Step(x,y)
end
