export add_type_D!, remove_type_D!

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
            if (norm(kink.i.vec-kink.k.vec) <= ex_radius)
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
            if (norm(kink.i.vec-kink.k.vec) <= ex_radius)
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
This function will find possibilites for the choice of the first Orbital that should be part of both Kinks after an add_type_D!-update,
when the two old orbitals of the old Kink are already selected.
To check possibility this will in particular check conservation laws and the occupation of the resulting fourth orb.
"""
function possible_new_orb1_D(occs, excite_orb1, excite_orb2, old_orb1, old_orb2)
    opportunities = filter(new_orb_1 -> in(find_fourth_orb_for_kink(new_orb_1, old_orb1, old_orb2),occs) && (new_orb_1 != find_fourth_orb_for_kink(new_orb_1, old_orb1, old_orb2)),
                    intersect!(sphere_with_same_spin(excite_orb1, dk = ex_radius), occs))

    return setdiff!(opportunities, Set([excite_orb1, excite_orb2, old_orb1, old_orb2]))
end

function add_type_D!(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64,Step}
    prop_prob = 1.0
    if isempty(c.kinks)
        return 1.0, Step()
    end

    old_τ, old_kink = rand(c.kinks)
    prop_prob *= 1.0/length(c.kinks)
    occs = occupations_at(c, old_τ)
    prop_prob *= 0.5 #left or right

    direction = rand([:left, :right])

    if direction == :left
        old_orbs = annihilators(old_kink)
        retained_orbs = creators(old_kink)
    else
        old_orbs = creators(old_kink)
        retained_orbs = annihilators(old_kink)
    end

    opportunities_new_orb1 = possible_new_orb1_D(occs, retained_orbs..., old_orbs...)

    if isempty(opportunities_new_orb1)
        return 1.0, Step()
    end

    new_orb1 = rand(opportunities_new_orb1)
    prop_prob *= 1.0/length(opportunities_new_orb1)
    new_orbs = (new_orb1, find_fourth_orb_for_kink(new_orb1, retained_orbs...))

    @assert(in(new_orbs[2], occs) & (new_orbs[1] != new_orbs[2]))

    if direction == :left
        τ_interval = old_τ - τ_prev_affecting(c.kinks, (retained_orbs..., new_orbs...), old_τ)
    else
        τ_interval = τ_next_affecting(c.kinks, (retained_orbs..., new_orbs...), old_τ ) - old_τ
    end

    if τ_interval < 0
        τ_interval +=1
    end
    @assert τ_interval > 0

    delta_τ = ImgTime(rand()*τ_interval)

    if direction == :left
        τ_new_kink = old_τ - delta_τ
    else
        τ_new_kink = old_τ + delta_τ
    end

    if τ_new_kink < 0
        change_occupations = true
        τ_new_kink += 1
    elseif τ_new_kink >= 1
        change_occupations = true
        τ_new_kink -= 1
    else
        change_occupations = false
    end

    while haskey(c.kinks, τ_new_kink)
        delta_τ = ImgTime(rand()*τ_interval)
        if dir == :left
            τ_new_kink = old_τ - delta_τ
        else
            τ_new_kink = old_τ + delta_τ
        end
        if τ_new_kink < 0
            change_occupations = true
            τ_new_kink += 1
        elseif τ_new_kink > 1
            change_occupations = true
            τ_new_kink -= 1
        else
            change_occupations = false
        end
    end

    if direction == :left
        τ_left = τ_new_kink
        τ_right = old_τ
    else
        τ_left = old_τ
        τ_right = τ_new_kink
    end

    delta_τ = abs(delta_τ)

    prop_prob *= 1.0/float(τ_interval)

    if change_occupations
        orbs_drop = new_orbs
        orbs_add = retained_orbs
    else
        orbs_drop = nothing
        orbs_add = nothing
    end

    prop_prob *= 0.5# shuffle orbitals of the changed kink
    if direction == :left
        kinks_add = ( τ_left => T4(retained_orbs..., new_orbs...), τ_right => T4(random_shuffle(new_orbs...)..., old_orbs...) )
    else
        kinks_add = ( τ_left => T4(old_orbs..., random_shuffle(new_orbs...)...), τ_right => T4(new_orbs..., retained_orbs...) )
    end

    @assert is_type_D(last(first(kinks_add)), last(last(kinks_add)))

    kinks_drop = (old_τ => old_kink,)

    # MC Step generated by this update
    Δ = Step(orbs_drop, kinks_drop, orbs_add, kinks_add)
    delta_di = ΔWdiag_element(m, e, c, retained_orbs..., new_orbs..., τ_left, τ_right)
    apply_step!(c,Δ)

    # calculate weight differance
    dw_off_diag = abs( Woffdiag_element(m, e, last(first(kinks_add))) * Woffdiag_element(m, e, last(last(kinks_add))) / Woffdiag_element(m, e, old_kink) )
    dw = e.β * dw_off_diag * exp(-(e.β * delta_τ*(energy(m, retained_orbs[1]) + energy(m, retained_orbs[2]) - energy(m, new_orbs[1]) - energy(m, new_orbs[2])) + e.β * delta_di))

    if direction == :left
        inverse_prop_prob = 0.5 / length( right_type_D_removable_pairs( c.kinks ) )
    else
        inverse_prop_prob = 0.5 / length( left_type_D_removable_pairs( c.kinks ) )
    end
    @assert delta_τ > 0
    @assert !isinf((inverse_prop_prob/prop_prob)*dw)
    return (inverse_prop_prob/prop_prob)*dw, Δ
end


function remove_type_D!(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64,Step}
    prop_prob = 0.5
    direction = rand([:left, :right])

    if direction == :left
        #removed kink left of changed kink
        opportunities = right_type_D_removable_pairs(c.kinks)
    else
        #removed kink right of changed kink
        opportunities = left_type_D_removable_pairs(c.kinks)
    end

    if isempty(opportunities)
        return 1.0, Step()
    end

    i_removed, i_changed = rand(opportunities)

    (removed_kink_τ, removed_kink) = c.kinks[i_removed]
    (changed_kink_τ, changed_kink) = c.kinks[i_changed]
    if direction == :left
        τ_left = removed_kink_τ
        τ_right = changed_kink_τ
    else
        τ_left = changed_kink_τ
        τ_right = removed_kink_τ
    end

    prop_prob *= 1.0/length(opportunities)

    if direction == :left
        (retained_orb1, retained_orb2, removed_orb1, removed_orb2) = orbs(removed_kink)
        (_, _, changed_orb1, changed_orb2) = orbs(changed_kink)

        @assert removed_kink.i != removed_kink.k

        add_kinks = (changed_kink_τ => T4(retained_orb1, retained_orb2, changed_orb1, changed_orb2),)
    else
        (removed_orb1, removed_orb2, retained_orb1, retained_orb2) = orbs(removed_kink)
        (changed_orb1, changed_orb2, _, _) = orbs(changed_kink)

        add_kinks = (changed_kink_τ => T4( changed_orb1, changed_orb2, retained_orb1, retained_orb2),)
    end

    # see if c.occupations change
    if τ_right < τ_left
        drop_orbs = (retained_orb1, retained_orb2)
        add_orbs = (removed_orb1, removed_orb2)
    else
        drop_orbs = nothing
        add_orbs = nothing
    end

    drop_kinks = (changed_kink_τ => changed_kink, removed_kink_τ => removed_kink)

    # MC Step generated by this update
    Δ = Step(drop_orbs, drop_kinks, add_orbs, add_kinks)
    apply_step!(c,Δ)

    #calculate reverse_prop_prob
    occs = occupations_at(c, changed_kink_τ)

    opportunities_reverse_new_orb1 = possible_new_orb1_D(occs, retained_orb1, retained_orb2, changed_orb1, changed_orb2)

    if direction == :left
        τ_interval = changed_kink_τ - τ_prev_affecting( c.kinks, (removed_orb1, removed_orb2, retained_orb1, retained_orb2), changed_kink_τ )
    else
        τ_interval = τ_next_affecting( c.kinks, (removed_orb1, removed_orb2, retained_orb1, retained_orb2), changed_kink_τ ) - changed_kink_τ
    end

    @assert in(removed_orb1, opportunities_reverse_new_orb1)

    if τ_interval < 0
        τ_interval +=1
    end

    inverse_prop_prob = (0.5/(length(c.kinks))) * (1.0/length(opportunities_reverse_new_orb1)) * (1.0/float(τ_interval)) * (1/2)

    dw_off_diag = abs( Woffdiag_element(m, e, last(first(drop_kinks))) *
                       Woffdiag_element(m, e, last(last(drop_kinks))) /
                       Woffdiag_element(m, e, last(first(add_kinks)) ))

    delta_di = ΔWdiag_element(m, e, c, retained_orb1, retained_orb2, removed_orb1, removed_orb2, τ_left, τ_right)
    delta_τ = float(τ_right - τ_left)

    if delta_τ < 0
        delta_τ +=1
    end

    dw = (1.0/e.β)*(1.0/dw_off_diag) * exp(e.β * delta_τ * ( energy(m, retained_orb1) + energy(m, retained_orb2) - energy(m, removed_orb1) - energy(m, removed_orb2)) + e.β * delta_di)

    @assert delta_τ > 0
    @assert !isnan((inverse_prop_prob/prop_prob) * dw)
    return (inverse_prop_prob/prop_prob)*dw, Δ
end

function isuseful(c::Configuration, up::typeof(add_type_D!))
    if isempty(c.kinks)
        return false
    else
        return true
    end
end

function isuseful(c::Configuration, up::typeof(remove_type_D!))
    if length(c.kinks) < 3
        return false
    else
        return true
    end
end
