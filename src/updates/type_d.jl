export add_type_D, remove_type_D

"""Returns True if left_kink and right_kink are entangled in a Type-D way
This does not check wether the two kinks are neighbouring"""
function is_type_D(left_kink::T4, right_kink::T4)
    if !(Set([left_kink.i, left_kink.j]) == Set([right_kink.k, right_kink.l])) &
        (Set([left_kink.k, left_kink.l]) == Set([right_kink.i, right_kink.j]))
        return true
    else
        return false
    end
end

"""Return pairs with indices of 'neighbouring' kinks that are Type-D entangled AND removable.
'neighbouring' refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-D-entanglement is oriented
#to the left of the first τ."""
function left_type_D_removable_pairs(ck)
    pairs_left = Set{Tuple{Int,Int}}()
    for (i,(τ,kink)) in enumerate(ck)
        i_left = prev_affecting(ck, orbs(kink) ,τ)
        if is_type_D(last(ck[i_left]), kink)
            if norm(kink.i.vec - kink.k.vec) <= ex_radius
                if kink.i.spin == kink.k.spin
                    push!(pairs_left, (i, i_left))
                end
            end
        end
    end
    return pairs_left
end

"""Return pairs with indices of 'neighbouring' Kinks that are Type-D-Entangeld AND removable.
'neighbouring' refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-D-entanglement is oriented
#to the right of the first τ."""
function right_type_D_removable_pairs(ck)
    pairs_right = Set{Tuple{Int,Int}}()
    for (i,(τ,kink)) in enumerate(ck)
        i_right = next_affecting(ck, orbs(kink), τ)
        if is_type_D(kink, last(ck[i_right]))
            if norm(kink.i.vec-kink.k.vec) <= ex_radius
                if kink.i.spin == kink.k.spin
                    push!(pairs_right, (i, i_right))
                end
            end
        end
    end
    return pairs_right
end


"""
    possible_new_orb1_D(occs, excite_orb1, excite_orb2, old_orb1, old_orb2)
This function will find possibilites for the choice of the first Orbital that should be part of both Kinks after an add_type_D-update,
when the two old orbitals of the old Kink are already selected.
To check possibility this will in particular check conservation laws and the occupation of the resulting fourth orb.
"""
function possible_new_orb1_D(occs, excite_orb1, excite_orb2, old_orb1, old_orb2)
    opportunities = filter(new_orb_1 -> in(find_fourth_orb_for_kink(new_orb_1, old_orb1, old_orb2),occs) && (new_orb_1 != find_fourth_orb_for_kink(new_orb_1, old_orb1, old_orb2)),
                    intersect!(sphere_with_same_spin(excite_orb1, dk = ex_radius), occs))

    return setdiff!(opportunities, Set([excite_orb1, excite_orb2, old_orb1, old_orb2]))
end

function add_type_D(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64,Step}
    prop_prob = 1.0
    if isempty(c.kinks)
        return 1.0, Step()
    end
    old_τ, old_kink = rand(c.kinks)
    prop_prob *= 1.0/length(c.kinks)
    occs = occupations_at(c, old_τ)
    prop_prob *= 0.5 #left or right
    if rand() >= 0.5
        #add kink left
        opportunities_new_orb1 = possible_new_orb1_D(occs, orbs(old_kink)...)
        if isempty(opportunities_new_orb1)
            return 1.0, Step()
        end
        new_orb1 = rand(opportunities_new_orb1)
        prop_prob *= 1.0/length(opportunities_new_orb1)
        new_orb2 = find_fourth_orb_for_kink(new_orb1, old_kink.k, old_kink.l)

        @assert(in(new_orb2, occs) & (new_orb1 != new_orb2))

        τ_interval = old_τ - τ_prev_affecting( c.kinks, Set([old_kink.i, old_kink.j, new_orb1, new_orb2]), old_τ )

        if τ_interval < 0
            τ_interval +=1
        end
        τ_new_kink = old_τ - ImgTime(rand()*τ_interval)
        if τ_new_kink < 0
            τ_new_kink += 1
            delta_τ = float(old_τ - τ_new_kink + 1)
        else
            delta_τ = float(old_τ - τ_new_kink)
        end
        @assert τ_interval > 0
        #no 2 kinks at same τ
        while haskey(c.kinks, τ_new_kink)
            τ_new_kink = old_τ - ImgTime(rand()*τ_interval)
            if τ_new_kink < 0
                τ_new_kink += 1
                delta_τ = float(old_τ - τ_new_kink + 1)
            else
                delta_τ = float(old_τ - τ_new_kink)
            end
        end

        prop_prob *= 1.0/float(τ_interval)

        # see if c.occupations change
        if τ_new_kink > old_τ  #new kink was added left of old kink
            orbs_drop = (new_orb1, new_orb2)
            orbs_add = (old_kink.i, old_kink.j)
        else
            orbs_drop = nothing
            orbs_add = nothing
        end

        prop_prob *= 0.5# shuffle creators of the changed kink
        kinks_add =(
                    τ_new_kink => T4(old_kink.i, old_kink.j, new_orb1, new_orb2),
                    # shuffle creators
                    old_τ => T4( random_shuffle( new_orb2, new_orb1 )..., old_kink.k, old_kink.l )
                    )
        @assert is_type_D(last(first(kinks_add)), last(last(kinks_add)))

        kinks_drop = (old_τ => old_kink,)

        # MC Step generated by this update
        Δ = Step(orbs_drop, kinks_drop, orbs_add, kinks_add)

        # calculate weight differance
        delta_di = ΔWdiag_element(m, e, c, old_kink.i, old_kink.j, new_orb1, new_orb2, τ_new_kink, old_τ)
        dw_off_diag = abs(Woffdiag_element(m, e, old_kink.i, old_kink.j, new_orb1, new_orb2) * Woffdiag_element(m, e, last(last(kinks_add))) /
                        Woffdiag_element(m, e, old_kink))
        dw = e.β * dw_off_diag * exp(-(e.β * delta_τ*(energy(m, old_kink.i) + energy(m, old_kink.j) - energy(m, new_orb1) - energy(m, new_orb2)) + e.β * delta_di))

        inverse_prop_prob = 0.5 / length( right_type_D_removable_pairs( apply_step(c,Δ).kinks ) )
    else
        #add kink right
        opportunities_new_orb1 = possible_new_orb1_D(occs, old_kink.k, old_kink.l, old_kink.i, old_kink.j)
        if isempty(opportunities_new_orb1)
            return 1.0, Step()
        end
        new_orb1 = rand(opportunities_new_orb1)
        prop_prob *= 1.0/length(opportunities_new_orb1)
        new_orb2 = find_fourth_orb_for_kink(new_orb1, old_kink.i, old_kink.j)

        @assert in(new_orb2, occs) & (new_orb1 != new_orb2)

        τ_interval = τ_next_affecting( c.kinks, Set([ old_kink.k, old_kink.l, new_orb1, new_orb2]), old_τ ) - old_τ
        if τ_interval < 0
            τ_interval +=1
        end
        τ_new_kink = old_τ + ImgTime(rand()*τ_interval)
        if τ_new_kink > 1
            τ_new_kink -= 1
            delta_τ = float(τ_new_kink - old_τ + 1)
        else
            delta_τ = float(τ_new_kink - old_τ)
        end

        #no 2 kinks at same τ
        @assert τ_interval > 0
        while haskey(c.kinks, τ_new_kink)
            τ_new_kink = old_τ + ImgTime(rand()*τ_interval)
            if τ_new_kink > 1
                τ_new_kink -= 1
                delta_τ = float(τ_new_kink - old_τ + 1)
            else
                delta_τ = float(τ_new_kink - old_τ)
            end
        end

        prop_prob *= 1.0/float(τ_interval)

        #see if c.occupations change
        if τ_new_kink < old_τ  #new kink was added right of old kink
            orbs_drop = (new_orb1, new_orb2)
            orbs_add = (old_kink.k, old_kink.l)
        else
            orbs_drop = nothing
            orbs_add = nothing
        end

        prop_prob *= 0.5# shuffle annihilators (!) of the changed kink
        kinks_add = (
                    τ_new_kink => T4(new_orb1, new_orb2, old_kink.k, old_kink.l),
                    # shuffle annihilators
                    old_τ => T4( old_kink.i, old_kink.j, random_shuffle(new_orb2, new_orb1 )... )
                    )
        @assert is_type_D(last(last(kinks_add)), last(first(kinks_add)))

        kinks_drop = (old_τ => old_kink,)
        # MC Step generated by this update
        Δ = Step(orbs_drop, kinks_drop, orbs_add, kinks_add)

        # calculate weight change
        delta_di = ΔWdiag_element(m, e, c, old_kink.k, old_kink.l, new_orb1, new_orb2, old_τ, τ_new_kink)
        dw_off_diag = abs( Woffdiag_element(m, e, last(first(kinks_add))) * Woffdiag_element(m, e, last(last(kinks_add))) / Woffdiag_element(m, e, old_kink) )

        dw = e.β * dw_off_diag * exp(-(e.β * delta_τ*(energy(m, old_kink.k) + energy(m, old_kink.l) - energy(m, new_orb1) - energy(m, new_orb2)) + e.β * delta_di))

        inverse_prop_prob = 0.5 / length( left_type_D_removable_pairs( apply_step(c,Δ).kinks ) )
    end

    @assert delta_τ > 0
    @assert !isinf((inverse_prop_prob/prop_prob)*dw)
    return (inverse_prop_prob/prop_prob)*dw, Δ
end

function remove_type_D(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64,Step}
    prop_prob = 0.5
    if rand() > 0.5
        #removed kink left of changed kink
        opportunities = right_type_D_removable_pairs(c.kinks)
        if isempty(opportunities)
            return 1.0, Step()
        end
 
        i_removed, i_changed = rand(opportunities)

        removed_kink_τ, removed_kink = c.kinks[i_removed]
        changed_kink_τ, changed_kink = c.kinks[i_changed]

        prop_prob *= 1.0/length(opportunities)
        #save those for later
        removed_orb1 = removed_kink.k
        removed_orb2 = removed_kink.l

        #change configuration
        @assert removed_orb1 != changed_kink.k

        #see if c.occupations change
        if removed_kink_τ > changed_kink_τ
            orbs_drop = (removed_kink.i, removed_kink.j)
            orbs_add = (removed_orb1, removed_orb2)
        else
            orbs_drop = nothing
            orbs_add = nothing
        end

        kinks_add = (changed_kink_τ => T4(removed_kink.i, removed_kink.j, changed_kink.k, changed_kink.l),)

        kinks_drop = (removed_kink_τ => removed_kink, changed_kink_τ => changed_kink, )

        # MC Step generated by this update
        Δ = Step(orbs_drop, kinks_drop, orbs_add, kinks_add)
        new_c = apply_step(c,Δ)

        # calculate reverse_prop_prob
        occs = occupations_at(new_c, changed_kink_τ)
        opportunities_reverse_new_orb1 = possible_new_orb1_D(occs, removed_kink.i, removed_kink.j, changed_kink.k, changed_kink.l)
        @assert in(removed_orb1, opportunities_reverse_new_orb1)
        τ_interval = changed_kink_τ - τ_prev_affecting( new_c.kinks, Set([removed_orb1, removed_orb2, last(first(kinks_add)).i, last(first(kinks_add)).j]), changed_kink_τ )

        if τ_interval < 0
            τ_interval +=1
        end

        # number of kinks is reduced by one with this update
        inverse_prop_prob = (0.5/(length(c.kinks)-1)) * (1.0/length(opportunities_reverse_new_orb1)) * (1.0/float(τ_interval)) * (1/2)

        # calculate weight change
        delta_di = ΔWdiag_element(m, e, new_c, last(first(kinks_add)).i, last(first(kinks_add)).j, removed_orb1, removed_orb2, removed_kink_τ, changed_kink_τ)

        dw_off_diag = abs( Woffdiag_element(m, e, removed_orb1, removed_orb2, last(first(kinks_add)).k, last(first(kinks_add)).l) * Woffdiag_element(m, e, last(first(kinks_add)).i, last(first(kinks_add)).j, removed_orb1, removed_orb2) / Woffdiag_element(m, e, last(first(kinks_add))) )

        delta_τ = float(changed_kink_τ - removed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end

        dw = (1.0/e.β)* (1.0/dw_off_diag) * exp(e.β * delta_τ * (energy(m, last(first(kinks_add)).i) + energy(m, last(first(kinks_add)).j) - energy(m,removed_orb1) - energy(m,removed_orb2)) + e.β * delta_di)

    else
        #removed kink right of changed kink
        opportunities = left_type_D_removable_pairs(c.kinks)
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

        # change configuration

        #see if c.occupations change
        if removed_kink_τ < changed_kink_τ# new kink was added right to old kink
            orbs_drop = (removed_kink.k, removed_kink.l)
            orbs_add = (removed_orb1, removed_orb2)
        else
            orbs_drop = nothing
            orbs_add = nothing
        end

        kinks_add = (changed_kink_τ => T4(changed_kink.i, changed_kink.j, removed_kink.k, removed_kink.l),)
        kinks_drop = (removed_kink_τ => removed_kink, changed_kink_τ => changed_kink)

        # MC Step generated by this update
        Δ = Step(orbs_drop, kinks_drop, orbs_add, kinks_add)
        new_c = apply_step(c,Δ)

        # calculate reverse_proposal probability
        occs = occupations_at(new_c, changed_kink_τ)
        opportunities_reverse_new_orb1 = possible_new_orb1_D(occs, last(first(kinks_add)).k, last(first(kinks_add)).l, last(first(kinks_add)).i, last(first(kinks_add)).j)
        @assert(in(removed_orb1, opportunities_reverse_new_orb1))
        τ_interval =  τ_next_affecting( new_c.kinks, Set([ last(first(kinks_add)).k, last(first(kinks_add)).l, removed_orb1, removed_orb2]), changed_kink_τ ) - changed_kink_τ
        if τ_interval < 0
            τ_interval +=1
        end

        # the number of kinks is reduced by one in this update
        inverse_prop_prob = (0.5/(length(c.kinks)-1)) * (1.0/length(opportunities_reverse_new_orb1)) *
                                 (1.0/float(τ_interval)) * (1/2)

        #calculate weight change
        delta_di = ΔWdiag_element(m, e, new_c, last(first(kinks_add)).k, last(first(kinks_add)).l, removed_orb1, removed_orb2, changed_kink_τ, removed_kink_τ)

        dw_off_diag = abs( Woffdiag_element(m, e, removed_orb1, removed_orb2,  last(first(kinks_add)).k, last(first(kinks_add)).l) *
                            Woffdiag_element(m, e, last(first(kinks_add)).i, last(first(kinks_add)).j, removed_orb1, removed_orb2) /
                            Woffdiag_element(m, e, last(first(kinks_add))) )

        delta_τ = float(removed_kink_τ - changed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1.0/e.β)*(1.0/dw_off_diag) * exp(e.β * delta_τ*(energy(m, last(first(kinks_add)).k) + energy(m, last(first(kinks_add)).l) - energy(m, removed_orb1) - energy(m, removed_orb2)) + e.β * delta_di)

    end

    @assert !isinf(dw)
    @assert delta_τ > 0
    @assert !isinf((inverse_prop_prob/prop_prob) * dw)

    return (inverse_prop_prob/prop_prob) * dw, Δ
end


function isuseful(c::Configuration, up::typeof(add_type_D))
    if isempty(c.kinks)
        return false
    else
        return true
    end
end

function isuseful(c::Configuration, up::typeof(remove_type_D))
    if length(c.kinks) < 3
        return false
    else
        return true
    end
end
