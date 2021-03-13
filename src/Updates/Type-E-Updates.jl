function add_type_E(c::Configuration, e::Ensemble) :: Tuple{Float64,Step}
    #After the updates the i and l komponents of both kinks will contain the old kink, wile the j and k components contain the old orbitals
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
        #now choose orbitals of the old kinks that shell also be part of the new kink
        prop_prob *= 0.5
        #choose exitation to of old kink to happen in first kink after update and exitation to happen second
        if rand() > 0.5
            new_kink_old_creator = last(old_kink).i
            changed_kink_old_creator = last(old_kink).j
        else
            new_kink_old_creator = last(old_kink).j
            changed_kink_old_creator = last(old_kink).i
        end
        prop_prob *= 0.5
        if rand() > 0.5
            new_kink_old_annihilator = last(old_kink).k
            changed_kink_old_annihilator = last(old_kink).l
        else
            new_kink_old_annihilator = last(old_kink).l
            changed_kink_old_annihilator = last(old_kink).k
        end
        #find occupied orb for creation of Type_E
        opportunities_new_kink_new_annihilator = intersect!(union!(get_sphere_with_same_spin(new_kink_old_creator, dk = ex_radius),
                                                                        get_sphere_with_same_spin(OrbitalHEG(new_kink_old_creator.vec, -new_kink_old_annihilator.spin), dk = ex_radius)), occs)
        delete!(opportunities_new_kink_new_annihilator, last(old_kink).i)
        delete!(opportunities_new_kink_new_annihilator, last(old_kink).j)
        if isempty(opportunities_new_kink_new_annihilator)
            return 1.0, Step()
        end
        new_kink_new_annihilator = rand(opportunities_new_kink_new_annihilator)
        prop_prob *= 1.0/length(opportunities_new_kink_new_annihilator)
        #calculate new creator
        if new_kink_new_annihilator.spin == new_kink_old_annihilator.spin
            new_kink_new_creator = OrbitalHEG(new_kink_old_annihilator.vec + (new_kink_new_annihilator.vec - new_kink_old_creator.vec), new_kink_old_creator.spin)
        else
            new_kink_new_creator = OrbitalHEG(new_kink_old_annihilator.vec + (new_kink_new_annihilator.vec - new_kink_old_creator.vec), -1 * new_kink_old_creator.spin)
        end
        if (in(new_kink_new_creator, occs) | (new_kink_new_creator == last(old_kink).k) | (new_kink_new_creator == last(old_kink).l))
            return 1.0, Step()
        end

        τ_Intervall = first(old_kink) - first(get_τ_borders(c, Set([new_kink_old_creator, new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator]), first(old_kink)))

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

        delta_di = get_change_diagonal_interaction(c, e, T4(new_kink_old_creator,new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator), τ_new_kink, first(old_kink))

        #change_Configuration
        #see if c.occupations change
        if τ_new_kink > first(old_kink)  #consider that new kink was added left of old kink
            @assert (!in(new_kink_new_creator, c.occupations))
            # change_occupations(c.occupations, T4(new_kink_old_creator,new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator))
            drop_orbs = Set([new_kink_new_annihilator, new_kink_old_annihilator])
            add_orbs = Set([new_kink_old_creator,new_kink_new_creator])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        #shuffle Indices TODO: is this necessary for ergodicy ?
        #shuffle changed kink
        prop_prob *= 1.0/16.0
        if rand() < 0.5
            if rand() < 0.5
                # shuffle both
                add_kink1 = first(old_kink) => T4(new_kink_new_annihilator, changed_kink_old_creator, changed_kink_old_annihilator, new_kink_new_creator)
            else
                # shuffle creators
                add_kink1 = first(old_kink) => T4(new_kink_new_annihilator, changed_kink_old_creator, new_kink_new_creator, changed_kink_old_annihilator)
            end
        else
            if rand() < 0.5
                # shuffle annihilators
                add_kink1 = first(old_kink) => T4(changed_kink_old_creator, new_kink_new_annihilator, changed_kink_old_annihilator, new_kink_new_creator)
            else
                # do not shuffle
                add_kink1 = first(old_kink) => T4(changed_kink_old_creator, new_kink_new_annihilator, new_kink_new_creator, changed_kink_old_annihilator)
            end
        end
        #shuffle new kink
        if rand() < 0.5
            if rand() < 0.5
                # shuffle both
                add_kink2 = τ_new_kink => T4(new_kink_new_creator, new_kink_old_creator, new_kink_old_annihilator, new_kink_new_annihilator)
            else
                # shuffle annihilator
                add_kink2 = τ_new_kink => T4(new_kink_old_creator,new_kink_new_creator, new_kink_old_annihilator, new_kink_new_annihilator)
            end
        else
            if rand() < 0.5
                # shuffle creators
                add_kink2 = τ_new_kink => T4(new_kink_new_creator, new_kink_old_creator, new_kink_new_annihilator, new_kink_old_annihilator)
            else
                # do not shuffle
                add_kink2 = τ_new_kink => T4(new_kink_old_creator, new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator)
            end
        end

        drop_kinks = old_kink

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks), Configuration(add_orbs, add_kink1, add_kink2))

        @assert(is_type_E(promote(c,Δ).kinks[τ_new_kink], promote(c,Δ).kinks[first(old_kink)]) != false)

        #calculate weight differance
        dw_off_diag = get_abs_offdiagonal_element(e,promote(c,Δ).kinks[τ_new_kink]) * get_abs_offdiagonal_element(e,promote(c,Δ).kinks[first(old_kink)]) /
                                                get_abs_offdiagonal_element(e,last(old_kink))
        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(get_energy(promote(c,Δ).kinks[τ_new_kink].i) + get_energy(promote(c,Δ).kinks[τ_new_kink].j)-
                                                         get_energy(promote(c,Δ).kinks[τ_new_kink].k) - get_energy(promote(c,Δ).kinks[τ_new_kink].l)) + delta_di))

        inverse_prop_prob = (1.0/length(get_right_type_E_removable_pairs(promote(c,Δ)))) * 0.5 * 0.25

        @assert(inverse_prop_prob != Inf)
    else
        #add kink right
        #now choose orbitals of the old kinks that shell also be part of the new kink
        prop_prob *= 0.5
        if rand() > 0.5
            new_kink_old_creator = last(old_kink).i
            changed_kink_old_creator = last(old_kink).j
        else
            new_kink_old_creator = last(old_kink).j
            changed_kink_old_creator = last(old_kink).i
        end
        prop_prob *= 0.5
        if rand() > 0.5
            new_kink_old_annihilator = last(old_kink).k
            changed_kink_old_annihilator = last(old_kink).l
        else
            new_kink_old_annihilator = last(old_kink).l
            changed_kink_old_annihilator = last(old_kink).k
        end
        #find occupied orb for creation of Type_E
        opportunities_new_kink_new_annihilator = setdiff!(union!(get_sphere_with_same_spin(new_kink_old_creator, dk = ex_radius),
                                                                        get_sphere_with_same_spin(OrbitalHEG(new_kink_old_creator.vec, -new_kink_old_annihilator.spin), dk = ex_radius)), occs)
        delete!(opportunities_new_kink_new_annihilator, last(old_kink).k)
        delete!(opportunities_new_kink_new_annihilator, last(old_kink).l)
        if isempty(opportunities_new_kink_new_annihilator)
            return 1.0, Step()
        end
        new_kink_new_annihilator = rand(opportunities_new_kink_new_annihilator)
        prop_prob *= 1.0/length(opportunities_new_kink_new_annihilator)

        if new_kink_new_annihilator.spin == new_kink_old_annihilator.spin
            new_kink_new_creator = OrbitalHEG(new_kink_old_annihilator.vec + (new_kink_new_annihilator.vec - new_kink_old_creator.vec), new_kink_old_creator.spin)
        else
            new_kink_new_creator = OrbitalHEG(new_kink_old_annihilator.vec + (new_kink_new_annihilator.vec - new_kink_old_creator.vec), -1 * new_kink_old_creator.spin)
        end

        if (!in(new_kink_new_creator, occs) | (new_kink_new_creator == last(old_kink).i) | (new_kink_new_creator == last(old_kink).j))
            return 1.0, Step()
        end

        τ_Intervall = last(get_τ_borders(c, Set([
                        new_kink_old_creator, new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator]),first(old_kink))) -
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
                                                          #Inverse new kink
        delta_di = get_change_diagonal_interaction(c, e, T4(new_kink_new_annihilator, new_kink_old_annihilator, new_kink_old_creator, new_kink_new_creator), first(old_kink), τ_new_kink)

        #change_Configuration
        #see if c.occupations change
        if τ_new_kink < first(old_kink)  #new kink was added right of old kink
            # change_occupations(c.occupations, T4(new_kink_new_annihilator, new_kink_old_annihilator, new_kink_old_creator, new_kink_new_creator))
            drop_orbs = Set([new_kink_old_creator, new_kink_new_creator])
            add_orbs = Set([new_kink_new_annihilator, new_kink_old_annihilator])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        # shuffle Indices TODO: is this necessary for ergodicy ?
        # shuffle changed kink
        prop_prob *= 1/16
        if rand() < 0.5
            if rand() < 0.5
                # shuffle both
                add_kink1 = first(old_kink) => T4(new_kink_new_annihilator, changed_kink_old_creator, changed_kink_old_annihilator, new_kink_new_creator)
            else
                # shuffle annihilators
                add_kink1 = first(old_kink) => T4(changed_kink_old_creator, new_kink_new_annihilator, changed_kink_old_annihilator, new_kink_new_creator)
            end
        else
            if rand() < 0.5
                # shuffle creators
                add_kink1 = first(old_kink) => T4(new_kink_new_annihilator, changed_kink_old_creator, new_kink_new_creator, changed_kink_old_annihilator)
            else
                add_kink1 = first(old_kink) => T4(changed_kink_old_creator, new_kink_new_annihilator, new_kink_new_creator, changed_kink_old_annihilator)
            end
        end
        # shuffle new kink
        if rand() < 0.5
            if rand() < 0.5
                # shuffle both
                add_kink2 = τ_new_kink => T4(new_kink_new_creator, new_kink_old_creator, new_kink_old_annihilator, new_kink_new_annihilator)
            else
                # shuffle annihilators
                add_kink2 = τ_new_kink => T4(new_kink_old_creator, new_kink_new_creator, new_kink_old_annihilator, new_kink_new_annihilator)
            end
        else
            if rand() < 0.5
                # shuffle creators
                add_kink2 = τ_new_kink => T4(new_kink_new_creator, new_kink_old_creator, new_kink_new_annihilator, new_kink_old_annihilator)
            else
                # do not shuffle
                add_kink2 = τ_new_kink => T4(new_kink_old_creator, new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator)
            end
        end

        drop_kinks = old_kink

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks), Configuration(add_orbs, add_kink1, add_kink2))

        dw_off_diag = get_abs_offdiagonal_element(e,promote(c,Δ).kinks[τ_new_kink]) * get_abs_offdiagonal_element(e,promote(c,Δ).kinks[first(old_kink)]) /
                                                get_abs_offdiagonal_element(e,last(old_kink))

        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(get_energy(promote(c,Δ).kinks[τ_new_kink].k) + get_energy(promote(c,Δ).kinks[τ_new_kink].l) -
                                                         get_energy(promote(c,Δ).kinks[τ_new_kink].i) - get_energy(promote(c,Δ).kinks[τ_new_kink].j)) + delta_di))

        inverse_prop_prob = (1.0/length(get_left_type_E_removable_pairs(promote(c,Δ)))) * 0.5 * 0.25

        @assert(inverse_prop_prob != Inf)
    end

    occupations(promote(c,Δ).occupations, promote(c,Δ).kinks)

    @assert(delta_τ > 0 )
    @assert(dw != Inf)
    @assert(prop_prob != 0)
    @assert(inverse_prop_prob != Inf)
    @assert(inverse_prop_prob != 0)
    return ((inverse_prop_prob/prop_prob)*dw), Δ
end

function remove_type_E(c::Configuration, e::Ensemble) :: Tuple{Float64,Step}
    prop_prob = 0.5
    if rand() > 0.5
        #removed kink left of changed kink
        opportunities = get_right_type_E_removable_pairs(c)
        if isempty(opportunities)
            return 1.0, Step()
        end
        removed_kink_tuple, changed_kink_tuple = rand(opportunities)
        removed_kink_τ =  first(removed_kink_tuple)
        changed_kink_τ =  first(changed_kink_tuple)
        prop_prob *= 1.0/length(opportunities)
        #safe thoose for later
        removed_kink = last(removed_kink_tuple)
        changed_kink_old = last(changed_kink_tuple)
        # if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it
        @assert(dot(removed_kink.i.vec-removed_kink.k.vec,
                    removed_kink.i.vec-removed_kink.k.vec) <= (ex_radius^2))


        #shuffle Indices TODO: is this necessary for ergodicy ?
        prop_prob *= 0.25
        if rand() < 0.5
            if rand() < 0.5
                # shuffle all
                add_kinks = (changed_kink_τ => T4(changed_kink_old.i, removed_kink.i, changed_kink_old.l, removed_kink.l),)
            else
                # shuffle annihilators
                add_kinks = (changed_kink_τ => T4(removed_kink.i, changed_kink_old.i, changed_kink_old.l, removed_kink.l),)
            end
        else
            if rand() < 0.5
                # shuffle creators
                add_kinks = (changed_kink_τ => T4(changed_kink_old.i, removed_kink.i, removed_kink.l, changed_kink_old.l),)
            else
                # do note shuffle
                add_kinks = (changed_kink_τ => T4(removed_kink.i, changed_kink_old.i, removed_kink.l, changed_kink_old.l),)
            end
        end

        # see if c.occupations change
        if removed_kink_τ > changed_kink_τ
            drop_orbs = Set([removed_kink.i,removed_kink.j])
            add_orbs = Set([removed_kink.k,removed_kink.l])
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
        opportunities_occ_orb_E = intersect!(union!(get_sphere_with_same_spin(OrbitalHEG(removed_kink.i.vec, 1), dk = ex_radius), get_sphere_with_same_spin(OrbitalHEG(removed_kink.i.vec, -1), dk = ex_radius)), occs)
        delete!(opportunities_occ_orb_E, promote(c,Δ).kinks[changed_kink_τ].i)
        delete!(opportunities_occ_orb_E, promote(c,Δ).kinks[changed_kink_τ].j)
        @assert(in(removed_kink.k,opportunities_occ_orb_E))
        τ_Intervall = changed_kink_τ - first(get_τ_borders(promote(c,Δ), Set([removed_kink.k, removed_kink.l,removed_kink.i, removed_kink.j]),changed_kink_τ))


        if τ_Intervall < 0
            τ_Intervall +=1
        end

        inverse_prop_prob = (0.5/length(promote(c,Δ).kinks)) * (1.0/length(opportunities_occ_orb_E)) *
                                 (1.0/Float64(τ_Intervall)) * (1/4) * (1/16)

        #calculate weight change
        delta_di = get_change_diagonal_interaction(promote(c,Δ), e, removed_kink, removed_kink_τ, changed_kink_τ)

        dw_off_diag = get_abs_offdiagonal_element(e,removed_kink) *
                        get_abs_offdiagonal_element(e,changed_kink_old) /
                            get_abs_offdiagonal_element(e,promote(c,Δ).kinks[changed_kink_τ])

        delta_τ = Float64(changed_kink_τ - removed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1.0/e.β)* (1.0/dw_off_diag) * exp(e.β * delta_τ * (get_energy(removed_kink.i) + get_energy(removed_kink.j) -
                                                                    get_energy(removed_kink.k) - get_energy(removed_kink.l)) + delta_di)

    else
        #removed kink right of changed kink
        opportunities = get_left_type_E_removable_pairs(c)
        if isempty(opportunities)
            return 1.0, Step()
        end
        removed_kink_tuple, changed_kink_tuple = rand(opportunities)
        removed_kink_τ = first(removed_kink_tuple)
        changed_kink_τ = first(changed_kink_tuple)
        prop_prob *= 1.0/length(opportunities)
        removed_kink = last(removed_kink_tuple)
        changed_kink_old = last(changed_kink_tuple)
        # if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it
        @assert(dot(removed_kink.i.vec-removed_kink.k.vec,
                    removed_kink.i.vec-removed_kink.k.vec) <= (ex_radius^2))

        #shuffle Indices TODO: is this necessary for ergodicy ?
        prop_prob *= 0.25
        if rand() < 0.5
            if rand() < 0.5
                # shuffle all
                add_kinks = (changed_kink_τ => T4(changed_kink_old.i, removed_kink.i, changed_kink_old.l, removed_kink.l),)
            else
                # shuffle annihilators
                add_kinks = (changed_kink_τ => T4(removed_kink.i, changed_kink_old.i, changed_kink_old.l, removed_kink.l),)
            end
        else
            if rand() < 0.5
                # shuffle creators
                add_kinks = (changed_kink_τ => T4(changed_kink_old.i, removed_kink.i, removed_kink.l, changed_kink_old.l),)
            else
                # do note shuffle
                add_kinks = (changed_kink_τ => T4(removed_kink.i, changed_kink_old.i, removed_kink.l, changed_kink_old.l),)
            end
        end

        #see if c.occupations change
        if removed_kink_τ < changed_kink_τ  #new kink was added right of old kink
            drop_orbs = Set([removed_kink.k,removed_kink.l])
            add_orbs = Set([removed_kink.i, removed_kink.j])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        drop_kinks = (changed_kink_τ => c.kinks[changed_kink_τ], removed_kink_τ => c.kinks[removed_kink_τ])

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))

        #calculate reverse_prop_prob
        occs = occupations(promote(c,Δ), changed_kink_τ)
        opportunities_unocc_orb_E = setdiff!(union!(get_sphere_with_same_spin(OrbitalHEG(removed_kink.i.vec, 1), dk = ex_radius), get_sphere_with_same_spin(OrbitalHEG(removed_kink.i.vec, -1), dk = ex_radius)), occs)
        delete!(opportunities_unocc_orb_E, promote(c,Δ).kinks[changed_kink_τ].k)
        delete!(opportunities_unocc_orb_E, promote(c,Δ).kinks[changed_kink_τ].l)
        @assert(in(removed_kink.k,opportunities_unocc_orb_E))
        τ_Intervall =  last(get_τ_borders(promote(c,Δ), Set([removed_kink.i, removed_kink.j,removed_kink.k, removed_kink.l]),changed_kink_τ)) - changed_kink_τ

        if τ_Intervall < 0
            τ_Intervall +=1
        end
        inverse_prop_prob = (0.5/length(promote(c,Δ).kinks)) * (1/length(opportunities_unocc_orb_E)) *
                                 (1.0/Float64(τ_Intervall)) * (1/4) * (1/16)

        #calculate weight change
        delta_di = get_change_diagonal_interaction(promote(c,Δ), e, T4(removed_kink.k,removed_kink.l, removed_kink.i,removed_kink.j), changed_kink_τ, removed_kink_τ)

        dw_off_diag = get_abs_offdiagonal_element(e,removed_kink) *
                        get_abs_offdiagonal_element(e,changed_kink_old) /
                            get_abs_offdiagonal_element(e,promote(c,Δ).kinks[changed_kink_τ])

        delta_τ = Float64(removed_kink_τ - changed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1.0/e.β)*(1.0/dw_off_diag) * exp(e.β * delta_τ*(get_energy(removed_kink.k) + get_energy(removed_kink.l) -
                                                            get_energy(removed_kink.i) - get_energy(removed_kink.j)) + delta_di)

    end

    occupations(promote(c,Δ).occupations, promote(c,Δ).kinks)


    @assert(dw != Inf)
    @assert(delta_τ > 0 )
    @assert(!isnan((inverse_prop_prob/prop_prob) * dw))

    return ((inverse_prop_prob/prop_prob) * dw), Δ
end
