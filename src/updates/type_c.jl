export add_type_C, remove_type_C

"""Returns True if left_kink and right_kink are entangled in a Type-C way.
This does not check wether the two kinks are neighbouring"""
function is_type_C(left_kink::T4, right_kink::T4)
    if (Set([left_kink.i, left_kink.j]) == Set([right_kink.k, right_kink.l])) &
        !(Set([left_kink.k, left_kink.l]) == Set([right_kink.i, right_kink.j]))
        return true
    else
        return false
    end
end

"""Return pairs with indices of 'neighbouring' kinks that are Type-C entangled AND removable.
'neighbouring' refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-C-entanglement is oriented
#to the left of the first τ."""
function left_type_C_removable_pairs(ck)
    pairs_left = Set{Tuple{Int,Int}}()
    for (i,(τ,kink)) in enumerate(ck)
        i_left = prev_affecting(ck, orbs(kink), τ)
        if is_type_C(last(ck[i_left]), kink)
            if norm(kink.i.vec - kink.k.vec) <= ex_radius
                if kink.i.spin == kink.k.spin
                    push!(pairs_left, (i, i_left))
                end
            end
        end
    end
    return pairs_left
end


"""Return pairs with indices of 'neighbouring' Kinks that are Type-C entangled AND removable.
'neighbouring' refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-C-entanglement is oriented
#to the right of the first τ."""
function right_type_C_removable_pairs(ck)
    pairs_right = Set{Tuple{Int,Int}}()
    for (i,(τ,kink)) in enumerate(ck)
        i_right = next_affecting(ck, orbs(kink), τ)
        if is_type_C(kink, last(ck[i_right]))
            if norm(kink.i.vec - kink.k.vec) <= ex_radius
                if kink.i.spin == kink.k.spin
                    push!(pairs_right, (i, i_right))
                end
            end
        end
    end
    return pairs_right
end


"""
    possible_new_orb1_C(occs, orb_c, orb_d)
This function will find possibilites for the choice of the first Orbital that should be part of both Kinks after an add_type_C-update,
when the two old orbitals of the old Kink are already selected.
To check possibility this will in particular check conservation laws and the occupation of the resulting fourth orb.
"""
function possible_new_orb1_C(occs, exite_orb1, exite_orb2, old_orb1, old_orb2)
    opprtunities = filter(new_orb_1 -> !in(find_fourth_orb_for_kink(new_orb_1, old_orb1, old_orb2),occs) && (new_orb_1 != find_fourth_orb_for_kink(new_orb_1, old_orb1, old_orb2)),
                    setdiff!(sphere_with_same_spin(exite_orb1, dk = ex_radius), occs))

    return setdiff!(opprtunities, Set([exite_orb1, exite_orb2, old_orb1, old_orb2]))
end

function add_type_C(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64,Step}
    prop_prob = 1.0
    if isempty(c.kinks)
        return 1.0, Step()
    end
    old_kink = rand(c.kinks)
    prop_prob *= 1.0/length(c.kinks)
    occs = occupations_at(c, first(old_kink))
    prop_prob *= 0.5 #left or right
    if rand() >= 0.5
        #add kink left
        opportunities_new_orb1 = possible_new_orb1_C(occs, last(old_kink).k, last(old_kink).l,last(old_kink).i, last(old_kink).j)
        if isempty(opportunities_new_orb1)
            return 1.0, Step()
        end
        new_orb1 = rand(opportunities_new_orb1)
        prop_prob *= 1.0/length(opportunities_new_orb1)
        new_orb2 = find_fourth_orb_for_kink(new_orb1, last(old_kink).i, last(old_kink).j)

        @assert(!in(new_orb2, occs) & (new_orb1 != new_orb2))

        τ_Intervall = first(old_kink) - τ_prev_affecting( c.kinks, Set([last(old_kink).k, last(old_kink).l, new_orb1, new_orb2]), first(old_kink) )

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

        # see if c.occupations change
        if τ_new_kink > first(old_kink)  #new kink was added left of old kink
            drop_orbs = Set([last(old_kink).k,last(old_kink).l])
            add_orbs = Set([new_orb1, new_orb2])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        prop_prob *= 0.5# shuffle annihilators of the changed kink
        add_kinks = (
                    τ_new_kink => T4(new_orb1, new_orb2, last(old_kink).k, last(old_kink).l),
                    # shuffle annihilators
                    first(old_kink) => T4( last(old_kink).i, last(old_kink).j, random_shuffle(new_orb2, new_orb1 )... )
                    )
        drop_kinks = (old_kink,)

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))

        # calculate weight differance
        delta_di = ΔWdiag_element(m, e, c, new_orb1, new_orb2, last(old_kink).k, last(old_kink).l, τ_new_kink, first(old_kink))
        dw_off_diag = abs( Woffdiag_element(m, e, new_orb1, new_orb2, last(old_kink).k, last(old_kink).l) * Woffdiag_element(m, e, last(last(add_kinks))) / Woffdiag_element(m, e, last(old_kink)) )
        dw = e.β * dw_off_diag * exp(-(e.β * delta_τ*(energy(m, new_orb1) + energy(m, new_orb2) - energy(m, last(old_kink).k) - energy(m, last(old_kink).l)) + e.β * delta_di))

        inverse_prop_prob = 0.5 / length( right_type_C_removable_pairs( apply_step(c,Δ).kinks ) )

        @assert !isinf(inverse_prop_prob) "change kink left: inverse_prop_prob = Inf"
    else
        # add kink right
        opportunities_new_orb1 = possible_new_orb1_C(occs, last(old_kink).i, last(old_kink).j, last(old_kink).k, last(old_kink).l)
        if isempty(opportunities_new_orb1)
            @assert(false)
            return 1.0, Step()
        end
        new_orb1 = rand(opportunities_new_orb1)
        prop_prob *= 1.0/length(opportunities_new_orb1)
        new_orb2 = find_fourth_orb_for_kink(new_orb1, last(old_kink).k, last(old_kink).l)

        @assert !in(new_orb2, occs) & (new_orb1 != new_orb2)
        τ_Intervall = τ_next_affecting( c.kinks , Set([last(old_kink).i, last(old_kink).j, new_orb1, new_orb2]), first(old_kink) ) - first(old_kink)

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

        # see if c.occupations change
        if τ_new_kink < first(old_kink)  #new kink was added right of old kink
            drop_orbs = Set([last(old_kink).i,last(old_kink).j])
            add_orbs = Set([new_orb1,new_orb2])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        prop_prob *= 0.5# shuffle creators of the changed kink
        add_kinks =(
                    τ_new_kink => T4(last(old_kink).i, last(old_kink).j, new_orb1, new_orb2),
                    # shuffle creators
                    first(old_kink) => T4( random_shuffle( new_orb2, new_orb1)... , last(old_kink).k, last(old_kink).l )
                    )
        drop_kinks = (old_kink,)

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))

        # weight factor
        delta_di = ΔWdiag_element(m, e, c, new_orb1, new_orb2, last(old_kink).i, last(old_kink).j, first(old_kink), τ_new_kink)
        dw_off_diag = abs( Woffdiag_element(m, e, last(old_kink).i, last(old_kink).j, new_orb1, new_orb2) * Woffdiag_element(m, e, last(last(add_kinks))) / Woffdiag_element(m, e, last(old_kink)) )
        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(energy(m, new_orb1) + energy(m, new_orb2) - energy(m, last(old_kink).i) - energy(m, last(old_kink).j)) + e.β * delta_di))

        inverse_prop_prob = 0.5 / length( left_type_C_removable_pairs( apply_step(c,Δ).kinks ) )

        @assert !isinf(inverse_prop_prob) "change kink right: inverse_prop_prob = Inf"
    end

    @assert !iszero(prop_prob)
    @assert !isinf(dw)

    @assert delta_τ > 0
    @assert !isnan((inverse_prop_prob/prop_prob)*dw)
    return (inverse_prop_prob/prop_prob)*dw, Δ
end

function remove_type_C(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64,Step}
    prop_prob = 0.5
    if rand() > 0.5
        #removed kink left of changed kink
        opportunities = right_type_C_removable_pairs(c.kinks)
        if isempty(opportunities)
            return 1.0, Step()
        end

        i_removed, i_changed = rand(opportunities)

        removed_kink_τ, removed_kink = c.kinks[i_removed]
        changed_kink_τ, changed_kink = c.kinks[i_changed]

        prop_prob *= 1.0/length(opportunities)
        #save those for later
        removed_orb1 = removed_kink.i
        removed_orb2 = removed_kink.j

        @assert removed_orb1 != removed_kink.k

        # see if c.occupations change
        if removed_kink_τ > changed_kink_τ
            drop_orbs = Set([removed_orb1,removed_orb2])
            add_orbs = Set([removed_kink.k,removed_kink.l])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        add_kinks = (changed_kink_τ => T4(changed_kink.i, changed_kink.j, removed_kink.k, removed_kink.l),)
        drop_kinks = (changed_kink_τ => changed_kink, removed_kink_τ => removed_kink)

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))
        new_c = apply_step(c,Δ)

        #calculate reverse_prop_prob
        occs = occupations_at(new_c, changed_kink_τ)
        opportunities_reverse_new_orb1 = possible_new_orb1_C(occs, last(first(add_kinks)).k, last(first(add_kinks)).l, last(first(add_kinks)).i, last(first(add_kinks)).j)# TODO: use directly the orbitals referenced via c.kinks[changed_kink_τ] ? for now this is left for clearity

        @assert(in(removed_orb1,opportunities_reverse_new_orb1))
        τ_Intervall = changed_kink_τ - τ_prev_affecting( new_c.kinks, Set([removed_orb1, removed_orb2, last(first(add_kinks)).k, last(first(add_kinks)).l]), changed_kink_τ )

        if τ_Intervall < 0
            τ_Intervall +=1
        end

        inverse_prop_prob = (0.5/(length(c.kinks)-1)) * (1.0/length(opportunities_reverse_new_orb1)) * (1.0/Float64(τ_Intervall)) * (1.0/2.0)# one kink is removed by this update

        #calculate weight change
        delta_di = ΔWdiag_element(m, e, new_c, removed_orb1, removed_orb2, removed_kink.k, removed_kink.l, removed_kink_τ, changed_kink_τ)

        dw_off_diag = abs( Woffdiag_element(m, e, removed_orb1, removed_orb2, last(first(add_kinks)).k, last(first(add_kinks)).l) *
                            Woffdiag_element(m, e, last(first(add_kinks)).i, last(first(add_kinks)).j, removed_orb1, removed_orb2) /
                            Woffdiag_element(m, e, last(first(add_kinks)) ))

        delta_τ = Float64(changed_kink_τ - removed_kink_τ)# TODO use float() ?
        if delta_τ < 0
            delta_τ +=1
        end

        dw = (1/e.β)* (1/dw_off_diag) * exp(e.β * delta_τ * (energy(m, removed_orb1) + energy(m, removed_orb2) -energy(m, last(first(add_kinks)).k) - energy(m, last(first(add_kinks)).l)) + e.β * delta_di)

    else
        # removed kink right of changed kink
        opportunities = left_type_C_removable_pairs(c.kinks)
        if isempty(opportunities)
            return 1.0, Step()
        end

        i_removed, i_changed = rand(opportunities)

        (removed_kink_τ, removed_kink) = c.kinks[i_removed]
        (changed_kink_τ, changed_kink) = c.kinks[i_changed]

        prop_prob *= 1.0/length(opportunities)
        #save thoose for later
        removed_orb1 = removed_kink.k
        removed_orb2 = removed_kink.l

        #change configuration
        add_kinks = (changed_kink_τ => T4(removed_kink.i, removed_kink.j, changed_kink.k, changed_kink.l),)

        #see if c.occupations change
        if removed_kink_τ < changed_kink_τ # new kink was added right of old kink
            drop_orbs = Set([removed_orb1,removed_orb2])
            add_orbs = Set([removed_kink.i, removed_kink.j])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        drop_kinks = (changed_kink_τ => changed_kink, removed_kink_τ => removed_kink)

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))
        new_c = apply_step(c,Δ)

        # calculate reverse_prop_prob
        occs = occupations_at(new_c, changed_kink_τ)
        opportunities_reverse_new_orb1 = possible_new_orb1_C(occs, last(first(add_kinks)).i, last(first(add_kinks)).j, last(first(add_kinks)).k, last(first(add_kinks)).l)

        @assert in(removed_orb1,opportunities_reverse_new_orb1)

        τ_Intervall =  τ_next_affecting( new_c.kinks, Set([last(first(add_kinks)).i, last(first(add_kinks)).j,removed_orb1, removed_orb2]), changed_kink_τ ) - changed_kink_τ

        if τ_Intervall < 0
            τ_Intervall +=1
        end
        inverse_prop_prob = (0.5/length(new_c.kinks)) * (1.0/length(opportunities_reverse_new_orb1)) * (1.0/Float64(τ_Intervall)) * (1/2)# TODO: use float() ?

        #calculate weight change
        delta_di = ΔWdiag_element(m, e, new_c, removed_orb1,removed_orb2, last(first(add_kinks)).i, last(first(add_kinks)).j, changed_kink_τ, removed_kink_τ)

        dw_off_diag = abs( Woffdiag_element(m, e, removed_orb1, removed_orb2,  last(first(add_kinks)).k, last(first(add_kinks)).l) *
                            Woffdiag_element(m, e, last(first(add_kinks)).i, last(first(add_kinks)).j, removed_orb1, removed_orb2) /
                            Woffdiag_element(m, e, last(first(add_kinks))) )

        delta_τ = Float64(removed_kink_τ - changed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1.0/e.β)*(1.0/dw_off_diag) * exp(e.β * delta_τ*(energy(m, removed_orb1) + energy(m, removed_orb2) - energy(m, last(first(add_kinks)).i) - energy(m, last(first(add_kinks)).j)) + e.β * delta_di)

    end

    @assert delta_τ > 0
    @assert !isnan((inverse_prop_prob/prop_prob) * dw)
    return (inverse_prop_prob/prop_prob)*dw, Δ
end

function isuseful(c::Configuration, up::typeof(add_type_C))
    if isempty(c.kinks)
        return false
    else
        return true
    end
end

function isuseful(c::Configuration, up::typeof(remove_type_C))
    if length(c.kinks) < 3
        return false
    else
        return true
    end
end
