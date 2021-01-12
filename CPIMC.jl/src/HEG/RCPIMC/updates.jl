const ex_radius = 2 #max Radius for exitation
function move_particle(c::Configuration, e::Ensemble) :: Tuple{Float64, Step}
    free_orbitals = get_non_interacting_orbs_of_set(c, c.occupations)
    if length(free_orbitals) == 0
        return 1.0, Step()
    else
        x = rand(get_non_interacting_orbs_of_set(c, c.occupations))
    end
    oe = get_non_interacting_orbs_of_set(c,setdiff!(get_sphere_with_same_spin(x, dk = ex_radius), c.occupations))

    #if there are no empty non interacting orbitals in neighbourhood make no change
    if length(oe) == 0
        return 1.0, Step()
    end
    y = rand(oe)
    @assert x != y "same Configuration proposed."

    delta_di = get_change_diagonal_interaction(c, e, T2(y,x), ImgTime(0), ImgTime(1))
    @assert delta_di != Inf
    # weight change
    dw = exp(-(e.β*(get_energy(y)-get_energy(x)) + delta_di))
    #println("kindiff:", e.β*(get_energy(y)-get_energy(x)))
    #println("wdiagdiff:", delta_di,"\n")

    # change Configuration
    Δ = Step(x, y)
    # delete!(c.occupations,x)
    # push!(c.occupations,y)

    # get orbitals for reverse update
    oe2 = get_non_interacting_orbs_of_set(promote(c,Δ),setdiff!(get_sphere_with_same_spin(y, dk = ex_radius), promote(c,Δ).occupations ))

    # quotient of proposal probabilities
    dv = length(oe)/length(oe2)
    @assert delta_di != Inf
    @assert (dv*dw) >= 0
    return dv*dw, Δ
end

function add_type_B(c::Configuration, e::Ensemble) :: Tuple{Float64, Step}
    #samplign propability
    prop_prob = 1
    #get first τ
    τ1 = ImgTime(rand())
    #do not allow two kinks at the same time
    while haskey(c.kinks, τ1)
        τ1 = ImgTime(rand())
    end
    occs = get_occupations_at(c, τ1)
    orb_c = rand(occs)
    prop_prob *= 1/e.N
    orb_d = rand(occs)
    while orb_d == orb_c
        orb_d = rand(occs)
    end
    prop_prob *= 1/(e.N-1)
    opportunities_orb_a = setdiff!(get_sphere_with_same_spin(orb_c, dk = ex_radius), occs)
    opportunities_orb_b = setdiff!(get_sphere_with_same_spin(orb_d, dk = ex_radius), occs)
    if (length(opportunities_orb_a) == 0) | (length(opportunities_orb_b) == 0)
        return 1.0, Step()
    end
    orb_a = rand(opportunities_orb_a)
    orb_b = OrbitalHEG((orb_c.vec-orb_a.vec) + orb_d.vec,orb_d.spin)
    @assert !((orb_a == orb_d) | (orb_b == orb_c))
    if (!in(orb_b,opportunities_orb_b) | (orb_a == orb_b))
        return 1.0, Step()
    end

    #We will change the proposal probability after we get τ2


    #get τ2, TODO: separate function
    borders = get_τ_borders(c, Set([orb_a,orb_b,orb_c,orb_d]),τ1)
    possible_τ2_interval = borders[2]-borders[1]
    if possible_τ2_interval < 0
        possible_τ2_interval += 1
    end
    τ2 = ImgTime(rand()*(possible_τ2_interval) + borders[1])
    if τ2 > 1
        τ2 -= 1
    end
    #do not allow two kinks at the same time
    while haskey(c.kinks, τ2)
        τ2 = ImgTime(rand()*(possible_τ2_interval) + borders[1])
        if τ2 > 1
            τ2 -= 1
        end
    end
    #See which of the two imiginary time point is the "left border of the intervall"
    #If there are no Kinks effecting any of the new Kinks orbitals then anyone
    #of the two τs can be firstτ
    if borders[2] == 1
        if rand() < 0.5
            firstτ = τ1
            lastτ = τ2
        else
            firstτ = τ2
            lastτ = τ1
        end
        prop_prob *= 0.5
    #If there are Kinks which τ is firstτ depends on which is the
    #left one inside the given Interval
    elseif τ1 > borders[1]
        if τ2 > borders[1]
            firstτ = min(τ1,τ2)
            lastτ = max(τ1,τ2)
        else
            firstτ = τ1
            lastτ = τ2
        end
    else
        if τ2 < borders[1]
            firstτ = min(τ1,τ2)
            lastτ = max(τ1,τ2)
        else
            firstτ = τ2
            lastτ = τ1
        end
    end
    delta_τ = lastτ - firstτ
    if delta_τ < 0
        delta_τ += 1
    end
    @assert(delta_τ > 0)
    @assert(delta_τ <= 1)


    #We do consider states that differ only threw the order off indices of kinks
    #as different states, that contribute all with the same weight with is already
    #blocked over all permutations (see function “get_abs_offdiagonal_element”),
    #therefore the updates where we end up with the same kinks but start building
    #the kink with a different Excitation will result in a different order off indices
    #and therefore considered a different Update (to compensate that we use a factor ¼
    #in the off diagonal Matrix element contribution) .

    #However we have to consider the different ways off getting to the same
    #update by choosing the imaginary times in a different order.

    #Therefore we modify the proposal_probability in the following way
    occs_τ2 = get_occupations_at(c, τ2)
    opportunities_orb_a_τ2 = setdiff!(get_sphere_with_same_spin(orb_c, dk = ex_radius), occs_τ2)
    @assert length(opportunities_orb_a_τ2) != 0
    prop_prob *= (1.0/length(opportunities_orb_a) + 1.0/length(opportunities_orb_a_τ2)) * 1.0/float(possible_τ2_interval)

    #calculate change in diagonal interaction energy from the current configuration c
    delta_di = get_change_diagonal_interaction(c, e, T4(orb_a,orb_b,orb_c,orb_d), firstτ, lastτ)

    #change configuration
    # c.kinks[firstτ] = T4(orb_a,orb_b,orb_c,orb_d)
    #shuffle index order of the second Kink
    if rand() < 0.5
        if rand() < 0.5
            #swap both
            # c.kinks[lastτ] = T4(orb_d,orb_c,orb_b,orb_a)
            add_kinks = (firstτ => T4(orb_a,orb_b,orb_c,orb_d),lastτ => T4(orb_d,orb_c,orb_b,orb_a))
        else
            #swap anihilators
            # c.kinks[lastτ] = T4(orb_c,orb_d,orb_b,orb_a)
            add_kinks = (firstτ => T4(orb_a,orb_b,orb_c,orb_d),lastτ => T4(orb_c,orb_d,orb_b,orb_a))
        end
    else
        if rand() < 0.5
            #swap creators
            # c.kinks[lastτ] = T4(orb_d,orb_c,orb_a,orb_b)
            add_kinks = (firstτ => T4(orb_a,orb_b,orb_c,orb_d), lastτ => T4(orb_d,orb_c,orb_a,orb_b))
        else
            #swap none
            # c.kinks[lastτ] = T4(orb_c,orb_d,orb_a,orb_b)
            add_kinks = (firstτ => T4(orb_a,orb_b,orb_c,orb_d), lastτ => T4(orb_c,orb_d,orb_a,orb_b))
        end
    end
    prop_prob *= 1/4

    #Look at wether the new pair of Kinks modifies the Occupations at τ = 0
    if firstτ > lastτ
        # change_occupations(c.occupations, T4(orb_a,orb_b,orb_c,orb_d))
        drop_orbs = Set([orb_c, orb_d])
        add_orbs = Set([orb_a, orb_b])
    else
        drop_orbs = nothing
        add_orbs = Set{basis(c)}()
    end

    # MC step generated by this update
    Δ = Step(drop_orbs, Configuration(add_orbs, add_kinks...))

    # quotient of proposal probabilities
    dv = (1/(length(promote(c,Δ).kinks)))/prop_prob# add 2 kinks proposed by this update

    # weight factor
    dw = ((e.β)^2) *
            get_abs_offdiagonal_element(e,promote(c,Δ),T4(orb_a,orb_b,orb_c,orb_d))^2 *
            exp(-((delta_τ)*e.β * (get_energy(orb_a) + get_energy(orb_b) -
                get_energy(orb_c) - get_energy(orb_d)) + delta_di))
    @assert (dv*dw) >= 0
    @assert isfinite(dv) "proposal probability in add_type_B! not finite."
    @assert isfinite(dw) "weight-change in add_type_B! is not finite."

    return dv*dw, Δ
end

function remove_type_B(c::Configuration, e::Ensemble) :: Tuple{Float64, Step}

    if length(c.kinks) == 0
        return 1.0, Step()
    end
    Kink1 = rand(c.kinks)
    # if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it
    if dot(last(Kink1).i.vec-last(Kink1).k.vec, last(Kink1).i.vec-last(Kink1).k.vec) > ex_radius^2
        return 1.0, Step()
    end
    prop_prob = 1/length(c.kinks)
    τ_Kink2 = last(get_τ_borders(c, Set([last(Kink1).i, last(Kink1).j, last(Kink1).k, last(Kink1).l]),first(Kink1)))
    Kink2 = τ_Kink2 => c.kinks[τ_Kink2]
    #look if Kinks are type-b-connected
    ijkl = Set([last(Kink1).i, last(Kink1).j, last(Kink1).k, last(Kink1).l])

    if ijkl == Set([last(Kink2).i, last(Kink2).j, last(Kink2).k, last(Kink2).l])
        # delete!(c.kinks, first(Kink1))
        # delete!(c.kinks, first(Kink2))
        #see if occupations at τ=0 are modified
        if first(Kink1) > first(Kink2)
            # change_occupations(c.occupations, last(Kink2))
            add_orbs = Set([last(Kink2).i, last(Kink2).j])
            drop_orbs = Set([last(Kink2).k, last(Kink2).l])
            delta_τ = first(Kink2)-first(Kink1) + 1
        else
            add_orbs = nothing
            drop_orbs = Set{basis(c)}()
            delta_τ = first(Kink2)-first(Kink1)
        end
        Δ = Step(Configuration(drop_orbs, Kink1, Kink2), add_orbs)
        @assert(delta_τ > 0)
        @assert(delta_τ <= 1)
        #calculate inverse prop_prob (see  add_type_B)
        borders = get_τ_borders(promote(c,Δ), ijkl ,first(Kink1))
        possible_τ2_interval = borders[2]-borders[1]
        if possible_τ2_interval < 0
            possible_τ2_interval = 1 + possible_τ2_interval
        end
        occs_τ_kink1 = get_occupations_at(promote(c,Δ), first(Kink1))
        occs_τ_kink2 = get_occupations_at(promote(c,Δ), first(Kink2))
        orb_a = last(Kink1).i
        orb_b = last(Kink1).j
        orb_c = last(Kink1).k
        orb_d = last(Kink1).l
        #See how prop_prob changes in the function add_type_B to understand this expression
        inverse_prop_prob = (1/e.N)*(1/(e.N-1)) *
            (1/length(setdiff!(get_sphere_with_same_spin(orb_c, dk = ex_radius), occs_τ_kink1))
                + 1/length(setdiff!(get_sphere_with_same_spin(orb_c, dk = ex_radius), occs_τ_kink2))) *
             1.0/float(possible_τ2_interval) * (1/4)
        if borders[2] == 1
            inverse_prop_prob *= 0.5
        end
        # quotient of proposal probabilities
        dv = inverse_prop_prob/prop_prob

        #calculate change in diagonal interaction energy
        delta_di = get_change_diagonal_interaction(promote(c, Δ), e, last(Kink1), first(Kink1), first(Kink2))

        # weight factor
        dw = (1.0/(e.β)^2) *
            (1.0/(get_abs_offdiagonal_element(e,promote(c,Δ),last(Kink1)))^2) *
                exp((delta_τ)*e.β * (get_energy(orb_a) +
                     get_energy(orb_b) - get_energy(orb_c) - get_energy(orb_d)) + delta_di)
        @assert (dv*dw) >= 0
        return dv*dw, Δ
    else
        return 1.0, Step()
    end
end


function change_type_B(c::Configuration, e::Ensemble) :: Tuple{Float64,Step}
    if length(c.kinks) == 0
        return 1.0, Step()
    end
    #print(length(c.kinks),"\n")
    Kink1 = rand(c.kinks)
    τ_Kink2 = last(get_τ_borders(c, Set([last(Kink1).i, last(Kink1).j]),first(Kink1)))
    Kink2 = τ_Kink2 => c.kinks[τ_Kink2]
    #look if Kinks are type-b-connected
    ijkl = Set([last(Kink1).i, last(Kink1).j, last(Kink1).k, last(Kink1).l])
    if ijkl != Set([last(Kink2).i, last(Kink2).j, last(Kink2).k, last(Kink2).l])
        return 1.0, Step()
    end
    occs = get_occupations_at(c, first(Kink1))

    opportunities = get_non_interacting_orbs_of_set_in_interval(
                        c,setdiff!(
                            get_sphere_with_same_spin(last(Kink1).i, dk = ex_radius
                            ), occs
                        ),first(Kink1),first(Kink2)
                    )
    delete!(opportunities, last(Kink1).k)
    delete!(opportunities, last(Kink1).l)
    if length(opportunities) == 0
        return 1.0, Step()
    end
    new_orb_i = rand(opportunities)
    new_orb_j = OrbitalHEG(last(Kink1).j.vec + last(Kink1).i.vec - new_orb_i.vec, last(Kink1).j.spin)
    if new_orb_i == new_orb_j
        return 1.0, Step()
    end
    if (!is_non_interacting_in_interval(c,new_orb_j,first(Kink1),first(Kink2)) |
        in(new_orb_j,occs))
        return 1.0, Step()
    else
        #calculate change in diagonal interaction energy from the old configuration c
        delta_di = get_change_diagonal_interaction(c, e, T4(new_orb_i, new_orb_j, last(Kink1).i, last(Kink1).j), first(Kink1), first(Kink2))

        #change occupations
        if first(Kink1) > first(Kink2)
            # change_occupations(c.occupations, T4(new_orb_i, new_orb_j, last(Kink1).i, last(Kink1).j))
            drop_orbs = Set([last(Kink1).i, last(Kink1).j])
            add_orbs = Set([new_orb_i, new_orb_j])
            delta_τ = first(Kink2)-first(Kink1) + 1
        else
            add_orbs = Set{basis(c)}()
            drop_orbs = Set{basis(c)}()
            delta_τ = first(Kink2)-first(Kink1)
        end


        dw = exp(-(e.β*delta_τ*(get_energy(new_orb_i) + get_energy(new_orb_j) -
                                    get_energy(last(Kink1).i) - get_energy(last(Kink1).j)) + delta_di)) *
            (get_abs_offdiagonal_element(e,promote(c, Step(drop_orbs,add_orbs)),T4(new_orb_i, new_orb_j, last(Kink1).k, last(Kink1).l))/
                        get_abs_offdiagonal_element(e,promote(c, Step(drop_orbs,add_orbs)),last(Kink1)))^2

        #change Kinks

        # c.kinks[first(Kink1)] = T4(new_orb_i, new_orb_j, last(Kink1).k, last(Kink1).l)
        drop_kinks = (Kink1,Kink2)
        if rand() <= 0.5
            # c.kinks[first(Kink2)] = T4(last(Kink2).i, last(Kink2).j, new_orb_i, new_orb_j,)
            add_kinks = (
                        first(Kink1) => T4(new_orb_i, new_orb_j, last(Kink1).k, last(Kink1).l),
                        first(Kink2) => T4(last(Kink2).i, last(Kink2).j, new_orb_i, new_orb_j)
                        )
        else
            # c.kinks[first(Kink2)] = T4(last(Kink2).i, last(Kink2).j, new_orb_j, new_orb_i,)
            add_kinks = (
                        first(Kink1) => T4(new_orb_i, new_orb_j, last(Kink1).k, last(Kink1).l),
                        first(Kink2) => T4(last(Kink2).i, last(Kink2).j, new_orb_j, new_orb_i,)
                        )
        end
        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))

        # calculate proposal probability
        change_occupations(occs, T4(new_orb_i, new_orb_j, last(Kink1).i, last(Kink1).j))
        opportunites_reverse = get_non_interacting_orbs_of_set_in_interval(
                                    promote(c,Δ),setdiff!(
                                        get_sphere_with_same_spin(new_orb_i, dk = ex_radius
                                        ), occs
                                    ),first(Kink1),first(Kink2)
                                )
        delete!(opportunites_reverse, last(Kink1).k)
        delete!(opportunites_reverse, last(Kink1).l)
        @assert (dw * length(opportunities)/length(opportunites_reverse)) >= 0
        return (dw * length(opportunities)/length(opportunites_reverse)),Δ
    end
end

function shuffle_indices(c::Configuration, e::Ensemble) :: Tuple{Float64, Step}
    if length(c.kinks) == 0
        return 1.0, Step()
    end
    kink = rand(c.kinks)
    if rand() > 0.5
        # c.kinks[first(kink)] = T4(last(kink).j,last(kink).i,last(kink).k,last(kink).l)
        return 1.0, Step(
                        Configuration(kink),
                        Configuration(first(kink) => T4(last(kink).j, last(kink).i, last(kink).k,last(kink).l))
                        )
    else
        # c.kinks[first(kink)] = T4(last(kink).i,last(kink).j,last(kink).l,last(kink).k)
        return 1.0, Step(
                        Configuration(kink),
                        Configuration(first(kink) => T4(last(kink).j, last(kink).i, last(kink).k,last(kink).l))
                        )
    end
end
