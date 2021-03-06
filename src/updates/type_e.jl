export add_type_E!, remove_type_E!

"""If left kink and right_kink are type-E-entangled it returns a tuple of the two kinks whos orbs
are sorted in a way that k an j of both kinks are the common orbitals.
Otherwise it returns false.
This does not check wether the two kinks are neighbouring"""
function is_type_E(left_kink::T4, right_kink::T4)
    c_orb1 = intersect(creators(left_kink), annihilators(right_kink))
    c_orb2 = intersect(annihilators(left_kink), creators(right_kink))
    if (length(c_orb1) == 1) & (length(c_orb2) == 1)
        c_orb1 = first(c_orb1)
        c_orb2 = first(c_orb2)
        noncommon_orb_leftk_left = first(setdiff(annihilators(left_kink), creators(right_kink)))
        noncommon_orb_rightk_right = first(setdiff(creators(right_kink), annihilators(left_kink)))
        noncommon_orb_leftk_right = first(setdiff(creators(left_kink), annihilators(right_kink)))
        noncommon_orb_rightk_left = first(setdiff(annihilators(right_kink), creators(left_kink)))
        @assert length(unique( (noncommon_orb_leftk_left, noncommon_orb_rightk_right, noncommon_orb_leftk_right, noncommon_orb_rightk_left ) ) ) == 4 " The given orbitals are not pairwise distinct. "
        @assert noncommon_orb_leftk_right.vec - noncommon_orb_leftk_left.vec == noncommon_orb_rightk_left.vec - noncommon_orb_rightk_right.vec
        return T4(noncommon_orb_leftk_right, c_orb1, c_orb2, noncommon_orb_leftk_left), T4(noncommon_orb_rightk_right, c_orb2, c_orb1, noncommon_orb_rightk_left)
    else
        return false
    end
end

"""Return a Set of tuples of pairs of two imaginary times and kinks of 'neighbouring' kinks that are Type-E entangled AND removable.
'neighbouring' refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-E-entanglement is oriented
#to the left of the first ??."""
function left_type_E_removable_pairs(ck)
    pairs_left = Set{Tuple{Pair{ImgTime,T4},Pair{ImgTime,T4}}}()
    for (i,(??,kink)) in enumerate(ck)
        kink_orbs = orbs(kink)
        i_left = prev_affecting(ck, kink_orbs, ??)
        sorted_kinks = is_type_E(last(ck[i_left]), kink)
        if sorted_kinks == false
            continue
        end
        i_right = next_affecting(ck, kink_orbs, ??)# TODO: this is only called to do the assertions
        @assert i_left != i_right
        @assert i_right != 0
        if norm(last(sorted_kinks).i.vec - last(sorted_kinks).k.vec) <= ex_radius
            push!(pairs_left, (?? => last(sorted_kinks), first(ck[i_left]) => first(sorted_kinks)))# TODO: how do we now that the first kink can be removed?
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
#to the right of the first ??."""
function right_type_E_removable_pairs(ck)
    pairs_right = Set{Tuple{Pair{ImgTime,T4},Pair{ImgTime,T4}}}()
    for (i,(??,kink)) in enumerate(ck)
        kink_orbs = orbs(kink)
        i_right = next_affecting(ck, kink_orbs, ??)

        sorted_kinks = is_type_E(kink, last(ck[i_right]))
        if sorted_kinks == false
            continue
        end
        i_left = prev_affecting(ck, kink_orbs, ??)# TODO: this is only called to do the assertions
        @assert i_left != i_right
        @assert i_right != 0
        if norm(first(sorted_kinks).i.vec - first(sorted_kinks).k.vec) <= ex_radius
            push!(pairs_right, (?? => first(sorted_kinks), first(ck[i_right]) => last(sorted_kinks)))
        end
    end
    return pairs_right
end

"""
    possible_new_kink_new_occ_orb(occs, new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator)
This function will find possibilites for the choice of the new Orbital that is occupied in occs and should be part of both Kinks after an add_type_E!-update,
when which orbital of the old-Kink is new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator and changed_kink_old_annihilator
is already selected.
To check possibility this will in particular check conservation laws and the occupation of the resulting fourth orb.
"""
function possible_new_kink_new_occ_orb(occs, new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator)
    possibilities = filter(new_kink_new_annihilator -> (!in(find_fourth_orb_for_kink(new_kink_old_creator, new_kink_old_annihilator, new_kink_new_annihilator),occs)
                                                        && (!in(find_fourth_orb_for_kink(new_kink_old_creator, new_kink_old_annihilator, new_kink_new_annihilator), (new_kink_old_annihilator, changed_kink_old_annihilator)))),
                    intersect!(union!(sphere_with_same_spin(new_kink_old_creator, dk = ex_radius),
                            sphere_with_same_spin(PlaneWave(new_kink_old_creator.vec, flip(new_kink_old_annihilator.spin)), dk = ex_radius)), occs))
    return setdiff!(possibilities, (new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator) )
end

"""
    possible_new_kink_new_occ_orb(occs, new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator)
This function will find possibilites for the choice of the new Orbital that is unoccupied in occs and should be part of both Kinks after an add_type_E!-update,
when which orbital of the old-Kink is new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator and changed_kink_old_annihilator
is already selected.
To check possibility this will in particular check conservation laws and the occupation of the resulting fourth orb.
"""
function possible_new_kink_new_unocc_orb(occs, new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator)
    possibilities = filter(new_kink_new_annihilator -> (in(find_fourth_orb_for_kink(new_kink_old_creator, new_kink_old_annihilator, new_kink_new_annihilator),occs)
                                                        && (!in(find_fourth_orb_for_kink(new_kink_old_creator, new_kink_old_annihilator, new_kink_new_annihilator), (new_kink_old_creator, changed_kink_old_creator)))),
                    setdiff!(union!(sphere_with_same_spin(new_kink_old_creator, dk = ex_radius),
                            sphere_with_same_spin(PlaneWave(new_kink_old_creator.vec, flip(new_kink_old_annihilator.spin)), dk = ex_radius)), occs))
    return setdiff!(possibilities, (new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator) )
end


function add_type_E!(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64,Step}
    #After the updates the i and l komponents of both kinks will contain the old kink, while the j and k components contain the old orbitals
    prop_prob = 1.0
    if isempty(c.kinks)
        return 1.0, Step()
    end
    old_kink_index = rand(1:length(c.kinks))# TODO: use this in order to pass the index of the removed kink to Step()
    old_kink = c.kinks[ old_kink_index ]
    prop_prob *= 1.0/length(c.kinks)
    occs = occupations_at(c, first(old_kink))
    prop_prob *= 0.5 # left or right
    direction = rand([:left, :right])
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
    if direction == :left
        opportunities_new_kink_new_annihilator = possible_new_kink_new_occ_orb(occs, new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator)
    else
        opportunities_new_kink_new_annihilator = possible_new_kink_new_unocc_orb(occs, new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator)
    end
    if isempty(opportunities_new_kink_new_annihilator)
        return 1.0, Step()
    end
    new_kink_new_annihilator = rand(opportunities_new_kink_new_annihilator)
    prop_prob *= 1.0/length(opportunities_new_kink_new_annihilator)
    #calculate new creator

    new_kink_new_creator = find_fourth_orb_for_kink(new_kink_old_creator, new_kink_old_annihilator, new_kink_new_annihilator)

    @assert(!in(new_kink_new_creator, orbs(last(old_kink))))

    if direction == :left
        @assert(!in(new_kink_new_creator, occs))
        ??_Intervall = first(old_kink) - ??_prev_affecting( c.kinks, (new_kink_old_creator, new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator), first(old_kink) )
    else
        @assert(in(new_kink_new_creator, occs))
        ??_Intervall = ??_next_affecting(c.kinks, (new_kink_old_creator, new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator), first(old_kink)) - first(old_kink)
    end

    if ??_Intervall < 0
        ??_Intervall +=1
    end
    @assert ??_Intervall > 0
    delta_?? = ImgTime(rand()*??_Intervall)
    if direction == :left
        ??_new_kink = first(old_kink) - delta_??
    else
        ??_new_kink = first(old_kink) + delta_??
    end
    if ??_new_kink > 1
        ??_new_kink -= 1
        change_orbs = true
    elseif ??_new_kink < 0
        ??_new_kink += 1
        change_orbs = true
    else
        change_orbs = false
    end
    #no 2 kinks at same ??
    while haskey(c.kinks, ??_new_kink)
        delta_?? = ImgTime(rand()*??_Intervall)
        if direction == :left
            ??_new_kink = first(old_kink) - delta_??
        else
            ??_new_kink = first(old_kink) + delta_??
        end
        if ??_new_kink > 1
            ??_new_kink -= 1
            change_orbs = true
        elseif ??_new_kink < 0
            ??_new_kink += 1
            change_orbs = true
        else
            change_orbs = false
        end
    end

    prop_prob *= 1.0/float(??_Intervall)

    # see if c.occupations change
    if change_orbs
        if direction == :left
            orbs_drop = (new_kink_new_annihilator, new_kink_old_annihilator)
            orbs_add = (new_kink_old_creator,new_kink_new_creator)
        else
            orbs_drop = (new_kink_old_creator, new_kink_new_creator)
            orbs_add = (new_kink_new_annihilator, new_kink_old_annihilator)
        end
    else
        orbs_drop = nothing
        orbs_add = nothing
    end

    # shuffle indices
    prop_prob *= 1.0/16.0
    # shuffle changed kink: shuffle creators and shuffle annihilators
    add_kink1 = first(old_kink) => T4( random_shuffle( new_kink_new_annihilator, changed_kink_old_creator )..., random_shuffle( changed_kink_old_annihilator, new_kink_new_creator )... )
    # shuffle new kink: shuffle creators and shuffle annihilators
    add_kink2 = ??_new_kink => T4( random_shuffle( new_kink_new_creator, new_kink_old_creator )..., random_shuffle( new_kink_old_annihilator, new_kink_new_annihilator )... )

    # MC Step generated by this update
    ?? = Step(orbs_drop, (old_kink,), orbs_add, (add_kink1, add_kink2))

    if direction == :left
        delta_di = ??Wdiag_element(m, e, c, new_kink_old_creator, new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator, ??_new_kink, first(old_kink))
    else
        delta_di = ??Wdiag_element(m, e, c, new_kink_new_annihilator, new_kink_old_annihilator, new_kink_old_creator, new_kink_new_creator, first(old_kink), ??_new_kink)
    end

    apply_step!(c,??)
    @assert (is_type_E(last(add_kink2), last(add_kink1)) != false)
    # calculate weight difference
    dw_off_diag = abs( Woffdiag_element(m, e, last(add_kink2)) * Woffdiag_element(m, e, last(add_kink1)) / Woffdiag_element(m, e, last(old_kink)) )
    if direction == :left
        dw = e.?? * dw_off_diag* exp(-(e.?? * delta_??*(energy(m, last(add_kink2).i) + energy(m, last(add_kink2).j)- energy(m, last(add_kink2).k) - energy(m, last(add_kink2).l)) + e.?? * delta_di))
        inverse_prop_prob = 0.125 / length( right_type_E_removable_pairs( c.kinks ) )# TODO: is there a faster way than to explicitly calculate all these pairs in order to get their number?
    else
        dw = e.?? * dw_off_diag* exp(-(e.?? * delta_??*(energy(m, last(add_kink2).k) + energy(m, last(add_kink2).l) - energy(m, last(add_kink2).i) - energy(m, last(add_kink2).j)) + e.?? * delta_di))
        inverse_prop_prob = 0.125 / length( left_type_E_removable_pairs( c.kinks ) )# TODO: is there a straighter way than to explicitly calculate all these pairs in order to get their number ?
    end
    @assert !isinf(inverse_prop_prob)
    @assert !iszero(prop_prob)
    @assert !iszero(inverse_prop_prob)
    @assert (  0 <= (inverse_prop_prob/prop_prob)*dw < Inf)
    return (inverse_prop_prob/prop_prob)*dw, ??
end

function remove_type_E!(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64,Step}
    prop_prob = 0.5
    direction = rand([:left, :right])
    #removed kink left of changed kink
    if direction == :left
        opportunities = right_type_E_removable_pairs(c.kinks)
    else
        opportunities = left_type_E_removable_pairs(c.kinks)
    end
    if isempty(opportunities)
        return 1.0, Step()
    end
    removed_kink_pair, changed_kink_pair = rand(opportunities)
    removed_kink_?? =  first(removed_kink_pair)
    changed_kink_?? =  first(changed_kink_pair)
    prop_prob *= 1.0/length(opportunities)

    # safe those for later
    removed_kink = last(removed_kink_pair)
    changed_kink_old = last(changed_kink_pair)
    # shuffle Indices
    prop_prob *= 0.25
    # shuffle creators and shuffle annihilators
    add_kink = changed_kink_?? => T4( random_shuffle( changed_kink_old.i, removed_kink.i )..., random_shuffle( changed_kink_old.l, removed_kink.l )...)
    if direction == :left
        ??_left = removed_kink_??
        ??_right = changed_kink_??
    else
        ??_left = changed_kink_??
        ??_right = removed_kink_??
    end
    # see if c.occupations change
    if ??_left > ??_right
        if direction == :left
            orbs_drop = (removed_kink.i,removed_kink.j)
            orbs_add = (removed_kink.k,removed_kink.l)
        else
            orbs_drop = (removed_kink.k,removed_kink.l)
            orbs_add = (removed_kink.i, removed_kink.j)
        end
    else
        orbs_drop = nothing
        orbs_add = nothing
    end

    kinks_drop = (changed_kink_?? => c.kinks[changed_kink_??], removed_kink_?? => c.kinks[removed_kink_??])

    # MC Step generated by this update
    ?? = Step(orbs_drop, kinks_drop, orbs_add, (add_kink,))

    apply_step!(c,??) # get proposed (new) configuration for inverse proposal probability

    # calculate reverse_prop_prob
    occs = occupations_at(c, changed_kink_??)
    if direction == :left
        opportunities_new_kink_new_annihilator = possible_new_kink_new_occ_orb(occs, removed_kink.i, changed_kink_old.i, removed_kink.l, changed_kink_old.l)
        @assert(in(removed_kink.k,opportunities_new_kink_new_annihilator))
        ??_Intervall = changed_kink_?? - ??_prev_affecting( c.kinks, (removed_kink.k, removed_kink.l,removed_kink.i, removed_kink.j), changed_kink_?? )
    else
        opportunities_new_kink_new_annihilator = possible_new_kink_new_unocc_orb(occs, removed_kink.i, changed_kink_old.i, removed_kink.l, changed_kink_old.l)
        @assert(in(removed_kink.k,opportunities_new_kink_new_annihilator))
        ??_Intervall =  ??_next_affecting( c.kinks, (removed_kink.i, removed_kink.j,removed_kink.k, removed_kink.l), changed_kink_?? ) - changed_kink_??
    end

    if ??_Intervall < 0
        ??_Intervall +=1
    end

    # the number of kinks is reduced by one with this update
    # TODO: explain all factors
    inverse_prop_prob = 0.5 / (length(c.kinks)) / length(opportunities_new_kink_new_annihilator) / float(??_Intervall) / 4 / 16

    # calculate weight change
    dw_off_diag = abs( Woffdiag_element(m, e, removed_kink) * Woffdiag_element(m, e, changed_kink_old) / Woffdiag_element(m, e, last(add_kink)) )

    delta_?? = float(??_right - ??_left)
    if delta_?? < 0
        delta_?? +=1
    end
    if direction == :left
        delta_di = ??Wdiag_element(m, e, c, removed_kink.i, removed_kink.j, removed_kink.k, removed_kink.l, removed_kink_??, changed_kink_??)
        dw = (1.0/e.??)* (1.0/dw_off_diag) * exp(e.?? * delta_?? * (energy(m,removed_kink.i) + energy(m,removed_kink.j) - energy(m,removed_kink.k) - energy(m,removed_kink.l)) + e.?? * delta_di)
    else
        delta_di = ??Wdiag_element(m, e, c, removed_kink.k, removed_kink.l, removed_kink.i,removed_kink.j, changed_kink_??, removed_kink_??)
        dw = (1.0/e.??)*(1.0/dw_off_diag) * exp(e.?? * delta_??*(energy(m, removed_kink.k) + energy(m, removed_kink.l) - energy(m, removed_kink.i) - energy(m, removed_kink.j)) + e.?? * delta_di)
    end
    @assert (  0 <= (inverse_prop_prob/prop_prob)*dw < Inf)
    return (inverse_prop_prob/prop_prob) * dw, ??
end

function isuseful(c::Configuration, up::typeof(add_type_E!))
    if isempty(c.kinks)
        return false
    else
        return true
    end
end

function isuseful(c::Configuration, up::typeof(remove_type_E!))
    if length(c.kinks) < 3
        return false
    else
        return true
    end
end
