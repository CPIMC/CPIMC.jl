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

    (old_kink_τ, old_kink) = rand(c.kinks)
    prop_prob *= 1.0/length(c.kinks)
    occs = occupations_at(c, old_kink_τ)
    prop_prob *= 0.5 #left or right

    direction = rand([:left, :right])
    
    if direction == :left
        #add kink left
        old_orbs = annihilators(old_kink)
        retained_orbs = creators(old_kink)
    else
        # add kink right
        old_orbs = creators(old_kink)
        retained_orbs = annihilators(old_kink)
    end
    
    opportunities_new_orb1 = possible_new_orb1_C(occs, old_orbs..., retained_orbs...)

    new_orbs = (x -> (x, find_fourth_orb_for_kink(x, retained_orbs...)))(rand(opportunities_new_orb1))
    
    if isempty(opportunities_new_orb1)
        return 1.0, Step()
    end
  
    prop_prob *= 1.0/length(opportunities_new_orb1)

    @assert(!in(new_orbs[2], occs) & (new_orbs[1] != new_orbs[2]))

    if direction == :left
        τ_interval = old_kink_τ - τ_prev_affecting(c.kinks, (old_orbs..., new_orbs...), old_kink_τ)
    else
        τ_interval = τ_next_affecting(c.kinks, (old_orbs..., new_orbs...), old_kink_τ) - old_kink_τ
    end

    if τ_interval < 0
        τ_interval +=1
    end
        
    if direction == :left
        τ_new_kink = old_kink_τ - ImgTime(rand()*τ_interval)
        if τ_new_kink < 0
            τ_new_kink += 1
            delta_τ = Float64(old_kink_τ - τ_new_kink + 1)
        else
            delta_τ = Float64(old_kink_τ - τ_new_kink)
        end
        @assert τ_interval > 0
        #no 2 kinks at same τ
        while haskey(c.kinks, τ_new_kink)
            τ_new_kink = old_kink_τ - ImgTime(rand()*τ_interval)
            if τ_new_kink < 0
                τ_new_kink += 1
                delta_τ = Float64(old_kink_τ - τ_new_kink + 1)
            else
                delta_τ = Float64(old_kink_τ - τ_new_kink)
            end
        end
    else
        τ_new_kink = old_kink_τ + ImgTime(rand()*τ_interval)
        if τ_new_kink > 1
            τ_new_kink -= 1
            delta_τ = Float64(τ_new_kink - old_kink_τ + 1)
        else
            delta_τ = Float64(τ_new_kink - old_kink_τ)
        end

        #no 2 kinks at same τ
        @assert τ_interval > 0
        while haskey(c.kinks, τ_new_kink)
            τ_new_kink = old_kink_τ + ImgTime(rand()*τ_interval)
            if τ_new_kink > 1
                τ_new_kink -= 1
                delta_τ = Float64(τ_new_kink - old_kink_τ + 1)
            else
                delta_τ = Float64(τ_new_kink - old_kink_τ)
            end
        end        
    end
        
    prop_prob *= 1.0/Float64(τ_interval)

    # see if c.occupations change
    if (direction == :left && τ_new_kink > old_kink_τ) || (direction == :right && τ_new_kink < old_kink_τ)
        drop_orbs = old_orbs
        add_orbs = new_orbs
    else
        drop_orbs = nothing
        add_orbs = nothing
    end

    prop_prob *= 0.5 # shuffle orbitals of the changed kink
    
    if direction == :left
        add_kinks = ( τ_new_kink => T4(new_orbs..., old_orbs...), old_kink_τ => T4(retained_orbs..., random_shuffle(new_orbs...)...) )
    else
        add_kinks = ( τ_new_kink => T4(old_orbs..., new_orbs...), old_kink_τ => T4(random_shuffle(new_orbs...)..., retained_orbs...) )
    end

    drop_kinks = (old_kink_τ => old_kink,)

    # MC Step generated by this update
    Δ = Step(drop_orbs, drop_kinks, add_orbs, add_kinks)

    # calculate weight differance
    if direction == :left
        delta_di = ΔWdiag_element(m, e, c, new_orbs..., old_orbs..., τ_new_kink, old_kink_τ)
        dw_off_diag = abs( Woffdiag_element(m, e, new_orbs..., old_orbs...) * Woffdiag_element(m, e, last(last(add_kinks))) / Woffdiag_element(m, e, old_kink) )
        inverse_prop_prob = 0.5 / length( right_type_C_removable_pairs( apply_step(c,Δ).kinks ) )
        @assert !isinf(inverse_prop_prob) "change kink left: inverse_prop_prob = Inf"
    else
        delta_di = ΔWdiag_element(m, e, c, new_orbs..., old_orbs..., old_kink_τ, τ_new_kink)
        dw_off_diag = abs( Woffdiag_element(m, e, old_orbs..., new_orbs...) * Woffdiag_element(m, e, last(last(add_kinks))) / Woffdiag_element(m, e, old_kink) )
        inverse_prop_prob = 0.5 / length( left_type_C_removable_pairs( apply_step(c,Δ).kinks ) )
        @assert !isinf(inverse_prop_prob) "change kink right: inverse_prop_prob = Inf"
    end
    
    dw = e.β * dw_off_diag * exp(-(e.β * delta_τ * (energy(m, new_orbs[1]) + energy(m, new_orbs[2]) - energy(m, old_orbs[1]) - energy(m, old_orbs[2])) + e.β * delta_di))

    @assert !iszero(prop_prob)
    @assert !isinf(dw)

    @assert delta_τ > 0
    @assert !isnan((inverse_prop_prob/prop_prob)*dw)
    return (inverse_prop_prob/prop_prob)*dw, Δ
end

function remove_type_C(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64,Step}
    prop_prob = 0.5
    direction = rand([:left, :right])
    
    if direction == :left
        #removed kink left of changed kink
        opportunities = right_type_C_removable_pairs(c.kinks)
    else
        #removed kink right of changed kink
        opportunities = left_type_C_removable_pairs(c.kinks)
    end
    
    if isempty(opportunities)
        return 1.0, Step()
    end
    
    i_removed, i_changed = rand(opportunities)

    (removed_kink_τ, removed_kink) = c.kinks[i_removed]
    (changed_kink_τ, changed_kink) = c.kinks[i_changed]
    
    prop_prob *= 1.0/length(opportunities)
    
    if direction == :left
        (removed_orb1, removed_orb2, retained_orb1, retained_orb2) = orbs(removed_kink)
        (changed_orb1, changed_orb2, _, _) = orbs(changed_kink)
        
        @assert removed_kink.i != removed_kink.k

        add_kinks = (changed_kink_τ => T4(changed_orb1, changed_orb2, retained_orb1, retained_orb2),)
    else
        (retained_orb1, retained_orb2, removed_orb1, removed_orb2) = orbs(removed_kink)
        (_, _, changed_orb1, changed_orb2) = orbs(changed_kink)

        add_kinks = (changed_kink_τ => T4(retained_orb1, retained_orb2, changed_orb1, changed_orb2),)
    end

    # see if c.occupations change
    if (direction == :left && removed_kink_τ > changed_kink_τ) || (direction == :right && removed_kink_τ < changed_kink_τ)
        drop_orbs = (removed_orb1, removed_orb2)
        add_orbs = (retained_orb1, retained_orb2)
    else
        drop_orbs = nothing
        add_orbs = nothing
    end
   
    drop_kinks = (changed_kink_τ => changed_kink, removed_kink_τ => removed_kink)

    # MC Step generated by this update
    Δ = Step(drop_orbs, drop_kinks, add_orbs, add_kinks)
    new_c = apply_step(c,Δ)

    #calculate reverse_prop_prob
    occs = occupations_at(new_c, changed_kink_τ)

    opportunities_reverse_new_orb1 = possible_new_orb1_C(occs, retained_orb1, retained_orb2, changed_orb1, changed_orb2)
   
    if direction == :left
        τ_interval = changed_kink_τ - τ_prev_affecting( new_c.kinks, (removed_orb1, removed_orb2, retained_orb1, retained_orb2), changed_kink_τ )
    else
        τ_interval = τ_next_affecting( new_c.kinks, (removed_orb1, removed_orb2, retained_orb1, retained_orb2), changed_kink_τ ) - changed_kink_τ
    end
    
    @assert in(removed_orb1, opportunities_reverse_new_orb1)

    if τ_interval < 0
        τ_interval +=1
    end

    inverse_prop_prob = (0.5/(length(c.kinks)-1)) * (1.0/length(opportunities_reverse_new_orb1)) * (1.0/float(τ_interval)) * (1/2)
    
    dw_off_diag = abs( Woffdiag_element(m, e, removed_orb1, removed_orb2, last(first(add_kinks)).k, last(first(add_kinks)).l) *
                       Woffdiag_element(m, e, last(first(add_kinks)).i, last(first(add_kinks)).j, removed_orb1, removed_orb2) /
                       Woffdiag_element(m, e, last(first(add_kinks)) ))

    if direction == :left
        delta_di = ΔWdiag_element(m, e, new_c, removed_orb1, removed_orb2, retained_orb1, retained_orb2, removed_kink_τ, changed_kink_τ)
        delta_τ = float(changed_kink_τ - removed_kink_τ)
    else
        delta_di = ΔWdiag_element(m, e, new_c, removed_orb1, removed_orb2, retained_orb1, retained_orb2, changed_kink_τ, removed_kink_τ)
        delta_τ = float(removed_kink_τ - changed_kink_τ)
    end
    
    if delta_τ < 0
        delta_τ +=1
    end
    
    dw = (1.0/e.β)*(1.0/dw_off_diag) * exp(e.β * delta_τ * (energy(m, removed_orb1) + energy(m, removed_orb2) - energy(m, retained_orb1) - energy(m, retained_orb2)) + e.β * delta_di)

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
