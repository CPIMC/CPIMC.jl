export add_type_E, remove_type_E

"""If left kink and right_kink are type-E-entangled it returns a tuple of the two kinks whos orbs
are sorted in a way that k an j of both kinks are the common orbitals.
Otherwise it returns false.
This does not check wether the two kinks are neighbouring"""
function is_type_E(left_kink::T4, right_kink::T4)
    c_orb1 = intersect!(Set([left_kink.i, left_kink.j]), Set([right_kink.k,right_kink.l]))
    c_orb2 = intersect!(Set([left_kink.k, left_kink.l]), Set([right_kink.i,right_kink.j]))
    if (length(c_orb1) == 1) & (length(c_orb2) == 1)
        c_orb1 = first(c_orb1)
        c_orb2 = first(c_orb2)
        noncommon_orb_leftk_left = first(setdiff!(Set([left_kink.k, left_kink.l]), Set([right_kink.i,right_kink.j])))
        noncommon_orb_rightk_right = first(setdiff!(Set([right_kink.i,right_kink.j]), Set([left_kink.k, left_kink.l])))
        noncommon_orb_leftk_right = first(setdiff!(Set([left_kink.i, left_kink.j]), Set([right_kink.k,right_kink.l])))
        noncommon_orb_rightk_left = first(setdiff!(Set([right_kink.k,right_kink.l]), Set([left_kink.i, left_kink.j])))
        @assert length(Set([noncommon_orb_leftk_left, noncommon_orb_rightk_right, noncommon_orb_leftk_right, noncommon_orb_rightk_left])) == 4 " The given orbitals are not pairwise distinct. "
        @assert noncommon_orb_leftk_right.vec - noncommon_orb_leftk_left.vec == noncommon_orb_rightk_left.vec - noncommon_orb_rightk_right.vec
        return T4(noncommon_orb_leftk_right, c_orb1, c_orb2, noncommon_orb_leftk_left), T4(noncommon_orb_rightk_right, c_orb2, c_orb1, noncommon_orb_rightk_left)
    else
        return false
    end
end

"""Return a Tuple of two imaginary times of 'neighbouring' kinks that are Type-E entangled AND removable.
'neighbouring' refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-E-entanglement is oriented
#to the left of the first τ."""
function left_type_E_removable_pairs(ck)
    pairs_left = Set{}()
    for (τ,kink) in ck
        kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
        τ_left,τ_right = τ_borders(ck, kink_orb_set ,τ)
        sorted_kinks = is_type_E(ck[τ_left], kink)
        if sorted_kinks == false
            continue
        end
        @assert τ_left != τ_right
        @assert τ_right != ImgTime(1)
        if norm(last(sorted_kinks).i.vec - last(sorted_kinks).k.vec) <= ex_radius
            push!(pairs_left, (τ => last(sorted_kinks),τ_left => first(sorted_kinks)))
        end
    end
    return pairs_left
end


"""Return a Tuple of 2 imaginary times of 'neighbouring' Kinks that are Type-E entangled AND removable.
'neighbouring' refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-E-entanglement is oriented
#to the right of the first τ."""
function right_type_E_removable_pairs(ck)
    pairs_right = Set{}()
    for (τ,kink) in ck
        kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
        τ_left,τ_right = τ_borders(ck, kink_orb_set ,τ)

        sorted_kinks = is_type_E(kink, ck[τ_right])
        if sorted_kinks == false
            continue
        end
        @assert τ_left != τ_right
        @assert τ_right != ImgTime(1)
        if norm(first(sorted_kinks).i.vec - first(sorted_kinks).k.vec) <= ex_radius
            push!(pairs_right, (τ => first(sorted_kinks),τ_right => last(sorted_kinks)))
        end
    end
    return pairs_right
end

"""
    possible_new_kink_new_occ_orb(occs, new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator)
This function will find possibilites for the choice of the new Orbital that is occupied in occs and should be part of both Kinks after an add_type_E-update,
when which orbital of the old-Kink is new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator and changed_kink_old_annihilator
is already selected.
To check possibility this will in particular check conservation laws and the occupation of the resulting fourth orb.
"""
function possible_new_kink_new_occ_orb(occs, new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator)
    possibilities = filter(new_kink_new_annihilator -> !in(find_fourth_orb_for_kink(new_kink_old_creator, new_kink_old_annihilator, new_kink_new_annihilator),occs),
                    intersect!(union!(sphere_with_same_spin(new_kink_old_creator, dk = ex_radius),
                            sphere_with_same_spin(PlaneWave(new_kink_old_creator.vec, flip(new_kink_old_annihilator.spin)), dk = ex_radius)), occs))
    return setdiff!(possibilities, Set([new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator]))
end

"""
    possible_new_kink_new_occ_orb(occs, new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator)
This function will find possibilites for the choice of the new Orbital that is unoccupied in occs and should be part of both Kinks after an add_type_E-update,
when which orbital of the old-Kink is new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator and changed_kink_old_annihilator
is already selected.
To check possibility this will in particular check conservation laws and the occupation of the resulting fourth orb.
"""
function possible_new_kink_new_unocc_orb(occs, new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator)
    possibilities = filter(new_kink_new_annihilator -> in(find_fourth_orb_for_kink(new_kink_old_creator, new_kink_old_annihilator, new_kink_new_annihilator),occs),
                    setdiff!(union!(sphere_with_same_spin(new_kink_old_creator, dk = ex_radius),
                            sphere_with_same_spin(PlaneWave(new_kink_old_creator.vec, flip(new_kink_old_annihilator.spin)), dk = ex_radius)), occs))
    return setdiff!(possibilities, Set([new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator]))
end


function add_type_E(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64,Step}
    #After the updates the i and l komponents of both kinks will contain the old kink, while the j and k components contain the old orbitals
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
        opportunities_new_kink_new_annihilator = possible_new_kink_new_occ_orb(occs, new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator)
        if isempty(opportunities_new_kink_new_annihilator)
            return 1.0, Step()
        end
        new_kink_new_annihilator = rand(opportunities_new_kink_new_annihilator)
        prop_prob *= 1.0/length(opportunities_new_kink_new_annihilator)
        #calculate new creator

        new_kink_new_creator = find_fourth_orb_for_kink(new_kink_old_creator, new_kink_old_annihilator, new_kink_new_annihilator)
        if (in(new_kink_new_creator, occs) | (new_kink_new_creator == last(old_kink).k) | (new_kink_new_creator == last(old_kink).l))
            return 1.0, Step()
        end

        τ_Intervall = first(old_kink) - τ_prev_affecting( c.kinks, Set([new_kink_old_creator, new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator]), first(old_kink) )

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
                delta_τ = Float64(first(old_kink) - τ_new_kink + 1)# TODO: use float() ?
            else
                delta_τ = Float64(first(old_kink) - τ_new_kink)# TODO: use float() ?
            end
        end

        prop_prob *= 1.0/Float64(τ_Intervall)# TODO: use float() ?

        # see if c.occupations change
        if τ_new_kink > first(old_kink) # consider that new kink was added left of old kink
            @assert !in(new_kink_new_creator, c.occupations)
            drop_orbs = Set([new_kink_new_annihilator, new_kink_old_annihilator])
            add_orbs = Set([new_kink_old_creator,new_kink_new_creator])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        # shuffle indices
        prop_prob *= 1.0/16.0
        # shuffle changed kink: shuffle creators and shuffle annihilators
        add_kink1 = first(old_kink) => T4( random_shuffle( new_kink_new_annihilator, changed_kink_old_creator )..., random_shuffle( changed_kink_old_annihilator, new_kink_new_creator )... )
        # shuffle new kink: shuffle creators and shuffle annihilators
        add_kink2 = τ_new_kink => T4( random_shuffle( new_kink_new_creator, new_kink_old_creator )..., random_shuffle( new_kink_old_annihilator, new_kink_new_annihilator )... )

        drop_kinks = old_kink

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks), Configuration(add_orbs, add_kink1, add_kink2))

        @assert (is_type_E(last(add_kink2), last(add_kink1)) != false)

        # calculate weight difference
        delta_di = ΔWdiag_element(m, e, c, new_kink_old_creator, new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator, τ_new_kink, first(old_kink))
        dw_off_diag = abs( Woffdiag_element(m, e, last(add_kink2)) * Woffdiag_element(m, e, last(add_kink1)) / Woffdiag_element(m, e, last(old_kink)) )
        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(energy(m, last(add_kink2).i) + energy(m, last(add_kink2).j)- energy(m, last(add_kink2).k) - energy(m, last(add_kink2).l)) + e.β * delta_di))

        inverse_prop_prob = 0.125 / length( right_type_E_removable_pairs( apply_step(c,Δ).kinks ) )# TODO: is there a faster way than to explicitly calculate all these pairs in order to get their number?

        @assert !isinf(inverse_prop_prob)
    else
        # add kink right
        # now choose orbitals of the old kinks that shell also be part of the new kink
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
        opportunities_new_kink_new_annihilator = possible_new_kink_new_unocc_orb(occs, new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator)
        if isempty(opportunities_new_kink_new_annihilator)
            return 1.0, Step()
        end
        new_kink_new_annihilator = rand(opportunities_new_kink_new_annihilator)
        prop_prob *= 1.0/length(opportunities_new_kink_new_annihilator)

        new_kink_new_creator = find_fourth_orb_for_kink(new_kink_old_creator, new_kink_old_annihilator, new_kink_new_annihilator)

        if !in(new_kink_new_creator, occs) | (new_kink_new_creator == last(old_kink).i) | (new_kink_new_creator == last(old_kink).j)
            return 1.0, Step()
        end

        τ_Intervall = τ_next_affecting(c.kinks, Set([new_kink_old_creator, new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator]), first(old_kink)) - first(old_kink)
        if τ_Intervall < 0
            τ_Intervall +=1
        end
        τ_new_kink = first(old_kink) + ImgTime(rand()*τ_Intervall)
        if τ_new_kink > 1
            τ_new_kink -= 1
            delta_τ = Float64(τ_new_kink - first(old_kink) + 1)# TODO: use float() ?
        else
            delta_τ = Float64(τ_new_kink - first(old_kink))# TODO: use float() ?
        end

        #no 2 kinks at same τ
        @assert τ_Intervall > 0
        while haskey(c.kinks, τ_new_kink)
            τ_new_kink = first(old_kink) + ImgTime(rand()*τ_Intervall)
            if τ_new_kink > 1
                τ_new_kink -= 1
                delta_τ = Float64(τ_new_kink - first(old_kink) + 1)# TODO: use float() ?
            else
                delta_τ = Float64(τ_new_kink - first(old_kink))# TODO: use float() ?
            end
        end

        prop_prob *= 1.0/Float64(τ_Intervall)# TODO: use float() ?

        # see if c.occupations change
        if τ_new_kink < first(old_kink)  #new kink was added right of old kink
            drop_orbs = Set([new_kink_old_creator, new_kink_new_creator])
            add_orbs = Set([new_kink_new_annihilator, new_kink_old_annihilator])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        # shuffle indices
        prop_prob *= 1/16
        # shuffle changed kink: shuffle creators and shuffle annihilators
        add_kink1 = first(old_kink) => T4( random_shuffle( new_kink_new_annihilator, changed_kink_old_creator )..., random_shuffle( changed_kink_old_annihilator, new_kink_new_creator )... )
        # shuffle new kink: shuffle creators and shuffle annihilators
        add_kink2 = τ_new_kink => T4( random_shuffle( new_kink_new_creator, new_kink_old_creator )..., random_shuffle( new_kink_old_annihilator, new_kink_new_annihilator )... )
        drop_kinks = old_kink

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks), Configuration(add_orbs, add_kink1, add_kink2))

        # calculate weight change

        # Inverse new kink TODO: what does this mean ?
        delta_di = ΔWdiag_element(m, e, c, new_kink_new_annihilator, new_kink_old_annihilator, new_kink_old_creator, new_kink_new_creator, first(old_kink), τ_new_kink)

        dw_off_diag = abs( Woffdiag_element(m, e, last(add_kink2)) * Woffdiag_element(m, e, last(add_kink1)) / Woffdiag_element(m, e, last(old_kink)) )

        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(energy(m, last(add_kink2).k) + energy(m, last(add_kink2).l) - energy(m, last(add_kink2).i) - energy(m, last(add_kink2).j)) + e.β * delta_di))

        inverse_prop_prob = 0.125 / length( left_type_E_removable_pairs( apply_step(c,Δ).kinks ) )# TODO: is there a straighter way than to explicitly calculate all these pairs in order to get their number ?

        @assert(inverse_prop_prob != Inf)
    end

    @assert !iszero(prop_prob)
    @assert !iszero(inverse_prop_prob)
    @assert (  0 <= (inverse_prop_prob/prop_prob)*dw < Inf)
    return (inverse_prop_prob/prop_prob)*dw, Δ
end

function remove_type_E(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64,Step}
    prop_prob = 0.5
    if rand() > 0.5
        #removed kink left of changed kink
        opportunities = right_type_E_removable_pairs(c.kinks)
        if isempty(opportunities)
            return 1.0, Step()
        end
        removed_kink_tuple, changed_kink_tuple = rand(opportunities)
        removed_kink_τ =  first(removed_kink_tuple)
        changed_kink_τ =  first(changed_kink_tuple)
        prop_prob *= 1.0/length(opportunities)
        # safe those for later
        removed_kink = last(removed_kink_tuple)
        changed_kink_old = last(changed_kink_tuple)
        #shuffle Indices
        prop_prob *= 0.25
        # shuffle creators and shuffle annihilators
        add_kink = changed_kink_τ => T4( random_shuffle( changed_kink_old.i, removed_kink.i )..., random_shuffle( changed_kink_old.l, removed_kink.l )... )

        # see if c.occupations change
        if removed_kink_τ > changed_kink_τ
            drop_orbs = Set([removed_kink.i,removed_kink.j])
            add_orbs = Set([removed_kink.k,removed_kink.l])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        drop_kinks = (changed_kink_τ => c.kinks[changed_kink_τ], removed_kink_τ => c.kinks[removed_kink_τ])

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kink))

        # calculate reverse_prop_prob
        occs = occupations_at(apply_step(c,Δ), changed_kink_τ)
        opportunities_occ_orb_E = possible_new_kink_new_occ_orb(occs, removed_kink.i, changed_kink_old.i, removed_kink.l, changed_kink_old.l)
        @assert(in(removed_kink.k,opportunities_occ_orb_E))
        τ_Intervall = changed_kink_τ - τ_prev_affecting( apply_step(c,Δ).kinks, Set([removed_kink.k, removed_kink.l,removed_kink.i, removed_kink.j]), changed_kink_τ )


        if τ_Intervall < 0
            τ_Intervall +=1
        end

        # the number of kinks is reduced by one with this update
        inverse_prop_prob = (0.5/(length(c.kinks)-1)) * (1.0/length(opportunities_occ_orb_E)) *
                                 (1.0/Float64(τ_Intervall)) * (1/4) * (1/16)# TODO: use float() ?

        # calculate weight change
        delta_di = ΔWdiag_element(m, e, apply_step(c,Δ), removed_kink.i, removed_kink.j, removed_kink.k, removed_kink.l, removed_kink_τ, changed_kink_τ)

        dw_off_diag = abs( Woffdiag_element(m, e, removed_kink) * Woffdiag_element(m, e, changed_kink_old) / Woffdiag_element(m, e, last(add_kink)) )

        delta_τ = Float64(changed_kink_τ - removed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1.0/e.β)* (1.0/dw_off_diag) * exp(e.β * delta_τ * (energy(m,removed_kink.i) + energy(m,removed_kink.j) - energy(m,removed_kink.k) - energy(m,removed_kink.l)) + e.β * delta_di)

    else
        #removed kink right of changed kink
        opportunities = left_type_E_removable_pairs(c.kinks)
        if isempty(opportunities)
            return 1.0, Step()
        end
        removed_kink_tuple, changed_kink_tuple = rand(opportunities)
        removed_kink_τ = first(removed_kink_tuple)
        changed_kink_τ = first(changed_kink_tuple)
        prop_prob *= 1.0/length(opportunities)
        removed_kink = last(removed_kink_tuple)
        changed_kink_old = last(changed_kink_tuple)
        #shuffle Indices
        prop_prob *= 0.25
        # shuffle creators and shuffle annihilators
        add_kink = changed_kink_τ => T4( random_shuffle( changed_kink_old.i, removed_kink.i )..., random_shuffle( changed_kink_old.l, removed_kink.l )... )
        # see if c.occupations change
        if removed_kink_τ < changed_kink_τ  #new kink was added right of old kink
            drop_orbs = Set([removed_kink.k,removed_kink.l])
            add_orbs = Set([removed_kink.i, removed_kink.j])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        drop_kinks = (changed_kink_τ => c.kinks[changed_kink_τ], removed_kink_τ => c.kinks[removed_kink_τ])

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kink))

        #calculate reverse_prop_prob
        occs = occupations_at(apply_step(c,Δ), changed_kink_τ)
        opportunities_unocc_orb_E = possible_new_kink_new_unocc_orb(occs, removed_kink.i, changed_kink_old.i, removed_kink.l, changed_kink_old.l)
        @assert(in(removed_kink.k,opportunities_unocc_orb_E))
        τ_Intervall =  τ_next_affecting( apply_step(c,Δ).kinks, Set([removed_kink.i, removed_kink.j,removed_kink.k, removed_kink.l]), changed_kink_τ ) - changed_kink_τ

        if τ_Intervall < 0
            τ_Intervall +=1
        end

        # the number of kinks is reduced by one with this update
        inverse_prop_prob = (0.5/(length(c.kinks)-1)) * (1/length(opportunities_unocc_orb_E)) * (1.0/Float64(τ_Intervall)) * (1/4) * (1/16)

        # calculate weight change
        delta_di = ΔWdiag_element(m, e, apply_step(c,Δ), removed_kink.k, removed_kink.l, removed_kink.i,removed_kink.j, changed_kink_τ, removed_kink_τ)

        dw_off_diag = abs( Woffdiag_element(m, e, removed_kink) * Woffdiag_element(m, e, changed_kink_old) / Woffdiag_element(m, e, last(add_kink)) )

        delta_τ = Float64(removed_kink_τ - changed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1.0/e.β)*(1.0/dw_off_diag) * exp(e.β * delta_τ*(energy(m, removed_kink.k) + energy(m, removed_kink.l) - energy(m, removed_kink.i) - energy(m, removed_kink.j)) + e.β * delta_di)

    end

    @assert (  0 <= (inverse_prop_prob/prop_prob)*dw < Inf)
    return (inverse_prop_prob/prop_prob) * dw, Δ
end

function isuseful(c::Configuration, up::typeof(add_type_E))
    if isempty(c.kinks)
        return false
    else
        return true
    end
end

function isuseful(c::Configuration, up::typeof(remove_type_E))
    if length(c.kinks) < 3
        return false
    else
        return true
    end
end
