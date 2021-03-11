function move_particle(c::Configuration, e::Ensemble) :: Tuple{Float64, Step}
    free_orbitals = get_non_interacting_orbs_of_set(c, c.occupations)
    if isempty(free_orbitals)
        return 1.0, Step()
    else
        x = rand(get_non_interacting_orbs_of_set(c, c.occupations))
    end
    oe = get_non_interacting_orbs_of_set(c,setdiff!(get_sphere_with_same_spin(x, dk = ex_radius), c.occupations))

    #if there are no empty non interacting orbitals in neighbourhood make no change
    if isempty(oe)
        return 1.0, Step()
    end
    y = rand(oe)
    @assert x != y "same Configuration proposed."

    delta_di = get_change_diagonal_interaction(c, e, T2(y,x), ImgTime(0), ImgTime(1))
    @assert delta_di != Inf
    # weight change
    dw = exp(-(e.β*(get_energy(y)-get_energy(x)) + delta_di))

    # MC Step generated by this update
    Δ = Step(x, y)

    # get orbitals for reverse update
    oe2 = get_non_interacting_orbs_of_set(promote(c,Δ),setdiff!(get_sphere_with_same_spin(y, dk = ex_radius), promote(c,Δ).occupations ))

    occupations(promote(c,Δ).occupations, promote(c,Δ).kinks)

    # quotient of proposal probabilities
    dv = length(oe)/length(oe2)
    @assert (dv*dw) >= 0
    return dv*dw, Δ
end
