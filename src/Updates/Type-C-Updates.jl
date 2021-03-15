function add_type_C(c::Configuration, e::Ensemble) :: Tuple{Float64,Step}
    prop_prob = 1.0
    if isempty(c.kinks)
        return 1.0, Step()
    end
    old_kink = rand(c.kinks)
    prop_prob *= 1.0/length(c.kinks)
    occs = occupations(c, first(old_kink))
    prop_prob *= 0.5 #left or right
    if rand() >= 0.5
        #add kink left
        opportunities_new_orb1 = setdiff!(get_sphere_with_same_spin(last(old_kink).k, dk = ex_radius), occs)
        delete!(opportunities_new_orb1, last(old_kink).k)
        delete!(opportunities_new_orb1, last(old_kink).l)
        if isempty(opportunities_new_orb1)
            return 1.0, Step()
        end
        new_orb1 = rand(opportunities_new_orb1)
        prop_prob *= 1.0/length(opportunities_new_orb1)
        new_orb2 = OrbitalHEG(last(old_kink).j.vec + (last(old_kink).i.vec - new_orb1.vec), last(old_kink).l.spin)
        if in(new_orb2, occs) | (new_orb1 == new_orb2)
            return 1.0, Step()
        end
        τ_Intervall = first(old_kink) - first(get_τ_borders(c, Set([last(old_kink).k, last(old_kink).l, new_orb1, new_orb2]), first(old_kink)))

        if τ_Intervall < 0
            τ_Intervall +=1
        end
        τ_new_kink = first(old_kink) - ImgTime(rand()*τ_Intervall)
        if τ_new_kink < 0
            τ_new_kink += 1
            delta_τ = Float64(first(old_kink) - τ_new_kink + 1)
        else
            delta_τ = Float64(first(old_kink) - τ_new_kink)
        end
        @assert τ_Intervall > 0
        #no 2 kinks at same τ
        while haskey(c.kinks, τ_new_kink)
            τ_new_kink = first(old_kink) - ImgTime(rand()*τ_Intervall)
            if τ_new_kink < 0
                τ_new_kink += 1
                delta_τ = Float64(first(old_kink) - τ_new_kink + 1)
            else
                delta_τ = Float64(first(old_kink) - τ_new_kink)
            end
        end

        prop_prob *= 1.0/Float64(τ_Intervall)

        delta_di = get_change_diagonal_interaction(c, e, T4(new_orb1,new_orb2,last(old_kink).k,last(old_kink).l), τ_new_kink, first(old_kink))

        #change_Configuration
        #see if c.occupations change
        if τ_new_kink > first(old_kink)  #new kink was added left of old kink
            # change_occupations(c.occupations, T4(new_orb1,new_orb2,last(old_kink).k,last(old_kink).l))
            drop_orbs = Set([last(old_kink).k,last(old_kink).l])
            add_orbs = Set([new_orb1, new_orb2])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        prop_prob *= 0.5# shuffle annihilators of the changed kink TODO: is this necessary for ergodicy ?
        if rand() < 0.5
            # shuffle
            add_kinks = (
                        τ_new_kink => T4(new_orb1, new_orb2, last(old_kink).k, last(old_kink).l),
                        first(old_kink) => T4(last(old_kink).i, last(old_kink).j, new_orb2, new_orb1)
                        )
        else
            # do not shuffle
            add_kinks = (
                        τ_new_kink => T4(new_orb1, new_orb2, last(old_kink).k, last(old_kink).l),
                        first(old_kink) => T4(last(old_kink).i, last(old_kink).j, new_orb1, new_orb2)
                        )
        end
        drop_kinks = (old_kink,)

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))

        #calculate weight differance
        dw_off_diag = get_abs_offdiagonal_element(e,promote(c,Δ).kinks[τ_new_kink]) * get_abs_offdiagonal_element(e,promote(c,Δ).kinks[first(old_kink)]) /
                                                get_abs_offdiagonal_element(e,last(old_kink))
        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(get_energy(new_orb1) + get_energy(new_orb2) -
                                    get_energy(last(old_kink).k) - get_energy(last(old_kink).l)) + delta_di))

        inverse_prop_prob = (1.0/length(get_right_type_C_pairs(promote(c,Δ)))) * 0.5

        @assert !isinf(inverse_prop_prob) "change kink left: inverse_prop_prob = Inf"

    else
        #add kink right
        opportunities_new_orb1 = setdiff!(get_sphere_with_same_spin(last(old_kink).i, dk = ex_radius), occs)
        delete!(opportunities_new_orb1, last(old_kink).k)
        delete!(opportunities_new_orb1, last(old_kink).l)
        if isempty(opportunities_new_orb1)
            return 1.0, Step()
        end
        new_orb1 = rand(opportunities_new_orb1)
        prop_prob *= 1.0/length(opportunities_new_orb1)
        new_orb2 = OrbitalHEG(last(old_kink).l.vec + (last(old_kink).k.vec - new_orb1.vec), last(old_kink).j.spin)
        if in(new_orb2, occs) | (new_orb1 == new_orb2)
            return 1.0, Step()
        end
        τ_Intervall = last(get_τ_borders(c, Set([
                        last(old_kink).i, last(old_kink).j, new_orb1, new_orb2]),first(old_kink))) -
                        first(old_kink)
        if τ_Intervall < 0
            τ_Intervall +=1
        end
        τ_new_kink = first(old_kink) + ImgTime(rand()*τ_Intervall)
        if τ_new_kink > 1
            τ_new_kink -= 1
            delta_τ = Float64(τ_new_kink - first(old_kink) + 1)
        else
            delta_τ = Float64(τ_new_kink - first(old_kink))
        end

        #no 2 kinks at same τ
        @assert τ_Intervall > 0
        while haskey(c.kinks, τ_new_kink)
            τ_new_kink = first(old_kink) + ImgTime(rand()*τ_Intervall)
            if τ_new_kink > 1
                τ_new_kink -= 1
                delta_τ = Float64(τ_new_kink - first(old_kink) + 1)
            else
                delta_τ = Float64(τ_new_kink - first(old_kink))
            end
        end

        prop_prob *= 1.0/Float64(τ_Intervall)

        delta_di = get_change_diagonal_interaction(c, e, T4(new_orb1,new_orb2,last(old_kink).i,last(old_kink).j), first(old_kink), τ_new_kink)

        #change_Configuration
        #see if c.occupations change
        if τ_new_kink < first(old_kink)  #new kink was added right of old kink
            # change_occupations(c.occupations, T4(new_orb1,new_orb2,last(old_kink).i,last(old_kink).j))
            drop_orbs = Set([last(old_kink).i,last(old_kink).j])
            add_orbs = Set([new_orb1,new_orb2])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end


        prop_prob *= 0.5# shuffle creators of the changed kink TODO: is this necessary for ergodicy ?
        if rand() < 0.5
            # shuffle
            add_kinks = (
                        τ_new_kink => T4(last(old_kink).i, last(old_kink).j, new_orb1, new_orb2),
                        first(old_kink) => T4(new_orb2, new_orb1, last(old_kink).k, last(old_kink).l)
                        )
        else
            # do not shuffle
            add_kinks = (
                        τ_new_kink => T4(last(old_kink).i, last(old_kink).j, new_orb1, new_orb2),
                        first(old_kink) => T4(new_orb1, new_orb2, last(old_kink).k, last(old_kink).l)
                        )
        end
        drop_kinks = (old_kink,)

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))



        dw_off_diag = get_abs_offdiagonal_element(e,promote(c,Δ).kinks[τ_new_kink]) * get_abs_offdiagonal_element(e,promote(c,Δ).kinks[first(old_kink)]) /
                                                get_abs_offdiagonal_element(e,last(old_kink))

        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(get_energy(new_orb1) + get_energy(new_orb2) -
                                    get_energy(last(old_kink).i) - get_energy(last(old_kink).j)) + delta_di))

        inverse_prop_prob = (1.0/length(get_left_type_C_pairs(promote(c,Δ)))) * 0.5

        @assert !isinf(inverse_prop_prob) "change kink right: inverse_prop_prob = Inf"
    end

    @assert(!iszero(prop_prob))
    @assert(!isinf(dw))

    @assert(delta_τ > 0 )
    @assert(!isnan((inverse_prop_prob/prop_prob)*dw))
    return ((inverse_prop_prob/prop_prob)*dw), Δ
end

function remove_type_C(c::Configuration, e::Ensemble) :: Tuple{Float64,Step}
    prop_prob = 0.5
    if rand() > 0.5
        #removed kink left of changed kink
        opportunities = get_right_type_C_pairs(c)
        if isempty(opportunities)
            return 1.0, Step()
        end
        removed_kink_τ, changed_kink_τ = rand(opportunities)
        prop_prob *= 1.0/length(opportunities)

        @assert(c.kinks[removed_kink_τ].i.spin == c.kinks[removed_kink_τ].k.spin)
        @assert ( dot(c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec,c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec) <= (ex_radius^2) ) "if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it"

        #safe thoose for later
        removed_orb1 = c.kinks[changed_kink_τ].k
        removed_orb2 = c.kinks[changed_kink_τ].l

        #change configuration
        # c.kinks[changed_kink_τ] = T4(c.kinks[changed_kink_τ].i, c.kinks[changed_kink_τ].j, c.kinks[removed_kink_τ].k,c.kinks[removed_kink_τ].l)
        # @assert removed_orb1 != c.kinks[changed_kink_τ].k
        @assert removed_orb1 != c.kinks[removed_kink_τ].k
        add_kinks = (changed_kink_τ => T4(c.kinks[changed_kink_τ].i, c.kinks[changed_kink_τ].j, c.kinks[removed_kink_τ].k,c.kinks[removed_kink_τ].l),)

        #see if c.occupations change
        if removed_kink_τ > changed_kink_τ
            # change_occupations(c.occupations, T4(c.kinks[changed_kink_τ].k,c.kinks[changed_kink_τ].l,removed_orb1,removed_orb2))
            drop_orbs = Set([removed_orb1,removed_orb2])
            add_orbs = Set([c.kinks[removed_kink_τ].k,c.kinks[removed_kink_τ].l])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        # delete!(c.kinks, removed_kink_τ)
        drop_kinks = (changed_kink_τ => c.kinks[changed_kink_τ], removed_kink_τ => c.kinks[removed_kink_τ])

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))

        #calculate reverse_prop_prob
        occs = occupations(promote(c,Δ), changed_kink_τ)
        opportunities_reverse_new_orb1 = setdiff!(get_sphere_with_same_spin(promote(c,Δ).kinks[changed_kink_τ].k, dk = ex_radius), occs)
        delete!(opportunities_reverse_new_orb1, promote(c,Δ).kinks[changed_kink_τ].k)
        delete!(opportunities_reverse_new_orb1, promote(c,Δ).kinks[changed_kink_τ].l)

        τ_Intervall = changed_kink_τ - first(get_τ_borders(promote(c,Δ), Set([removed_orb1, removed_orb2,
                         promote(c,Δ).kinks[changed_kink_τ].k, promote(c,Δ).kinks[changed_kink_τ].l]),changed_kink_τ))


        if τ_Intervall < 0
            τ_Intervall +=1
        end

        inverse_prop_prob = (0.5/length(promote(c,Δ).kinks)) * (1.0/length(opportunities_reverse_new_orb1)) *
                                 (1.0/Float64(τ_Intervall)) * (1.0/2.0)# TODO: (1.0/2.0) = 0.5

        #calculate weight change
        delta_di = get_change_diagonal_interaction(promote(c,Δ), e, T4(removed_orb1,removed_orb2, promote(c,Δ).kinks[changed_kink_τ].k, promote(c,Δ).kinks[changed_kink_τ].l), removed_kink_τ, changed_kink_τ)

        dw_off_diag = get_abs_offdiagonal_element(e,T4(removed_orb1, removed_orb2,  promote(c,Δ).kinks[changed_kink_τ].k, promote(c,Δ).kinks[changed_kink_τ].l)) *
                        get_abs_offdiagonal_element(e,T4(promote(c,Δ).kinks[changed_kink_τ].i, promote(c,Δ).kinks[changed_kink_τ].j, removed_orb1, removed_orb2)) /
                            get_abs_offdiagonal_element(e,promote(c,Δ).kinks[changed_kink_τ])

        delta_τ = Float64(changed_kink_τ - removed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1/e.β)* (1/dw_off_diag) * exp(e.β * delta_τ * (get_energy(removed_orb1) + get_energy(removed_orb2) -
                                    get_energy(promote(c,Δ).kinks[changed_kink_τ].k) - get_energy(promote(c,Δ).kinks[changed_kink_τ].l)) + delta_di)

    else
        #removed kink right of changed kink
        opportunities = get_left_type_C_pairs(c)
        if isempty(opportunities)
            return 1.0, Step()
        end
        removed_kink_τ, changed_kink_τ = rand(opportunities)
        prop_prob *= 1.0/length(opportunities)

        @assert(c.kinks[removed_kink_τ].i.spin == c.kinks[removed_kink_τ].k.spin)
        @assert (dot(c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec,c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec) <= (ex_radius^2)) "if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it"

        #safe thoose for later
        removed_orb1 = c.kinks[changed_kink_τ].i
        removed_orb2 = c.kinks[changed_kink_τ].j


        #change configuration
        # c.kinks[changed_kink_τ] = T4(c.kinks[removed_kink_τ].i, c.kinks[removed_kink_τ].j, c.kinks[changed_kink_τ].k,c.kinks[changed_kink_τ].l)
        add_kinks = (changed_kink_τ => T4(c.kinks[removed_kink_τ].i, c.kinks[removed_kink_τ].j, c.kinks[changed_kink_τ].k,c.kinks[changed_kink_τ].l),)

        #see if c.occupations change
        if removed_kink_τ < changed_kink_τ  #new kink was added right of old kink
            # change_occupations(c.occupations, T4(c.kinks[changed_kink_τ].i,c.kinks[changed_kink_τ].j,removed_orb1,removed_orb2))
            drop_orbs = Set([removed_orb1,removed_orb2])
            add_orbs = Set([c.kinks[removed_kink_τ].i, c.kinks[removed_kink_τ].j])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        # delete!(c.kinks, removed_kink_τ)
        drop_kinks = (changed_kink_τ => c.kinks[changed_kink_τ], removed_kink_τ => c.kinks[removed_kink_τ])

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))

        #calculate reverse_prop_prob
        occs = occupations(promote(c,Δ), changed_kink_τ)
        opportunities_reverse_new_orb1 = setdiff!(get_sphere_with_same_spin(promote(c,Δ).kinks[changed_kink_τ].i, dk = ex_radius), occs)
        delete!(opportunities_reverse_new_orb1, promote(c,Δ).kinks[changed_kink_τ].k)
        delete!(opportunities_reverse_new_orb1, promote(c,Δ).kinks[changed_kink_τ].l)

        τ_Intervall =  last(get_τ_borders(promote(c,Δ), Set([promote(c,Δ).kinks[changed_kink_τ].i, promote(c,Δ).kinks[changed_kink_τ].j,
                                                    removed_orb1, removed_orb2]),changed_kink_τ)) - changed_kink_τ

        if τ_Intervall < 0
            τ_Intervall +=1
        end
        inverse_prop_prob = (0.5/length(promote(c,Δ).kinks)) * (1.0/length(opportunities_reverse_new_orb1)) *
                                 (1.0/Float64(τ_Intervall)) * (1/2)# TODO: (1/2) = 0.5

        #calculate weight change
        delta_di = get_change_diagonal_interaction(promote(c,Δ), e, T4(removed_orb1,removed_orb2, promote(c,Δ).kinks[changed_kink_τ].i,promote(c,Δ).kinks[changed_kink_τ].j), changed_kink_τ, removed_kink_τ)

        dw_off_diag = get_abs_offdiagonal_element(e,T4(removed_orb1, removed_orb2,  promote(c,Δ).kinks[changed_kink_τ].k, promote(c,Δ).kinks[changed_kink_τ].l)) *
                        get_abs_offdiagonal_element(e,T4(promote(c,Δ).kinks[changed_kink_τ].i, promote(c,Δ).kinks[changed_kink_τ].j, removed_orb1, removed_orb2)) /
                            get_abs_offdiagonal_element(e,promote(c,Δ).kinks[changed_kink_τ])

        delta_τ = Float64(removed_kink_τ - changed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1.0/e.β)*(1.0/dw_off_diag) * exp(e.β * delta_τ*(get_energy(removed_orb1) + get_energy(removed_orb2) -
                                    get_energy(promote(c,Δ).kinks[changed_kink_τ].i) - get_energy(promote(c,Δ).kinks[changed_kink_τ].j)) + delta_di)

    end

    @assert(delta_τ > 0 )
    @assert(!isnan((inverse_prop_prob/prop_prob) * dw))
    return ((inverse_prop_prob/prop_prob) * dw), Δ
end
