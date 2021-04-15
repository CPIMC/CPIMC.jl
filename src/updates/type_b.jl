export add_type_B, remove_type_B, change_type_B

"""Returns True if left_kink and right_kink are entangled in a Type-B way.
This does not check wether the two kinks are neighbouring"""
function is_type_B(left_kink::T4, right_kink::T4)
  if (Set([left_kink.i, left_kink.j]) == Set([right_kink.k, right_kink.l])) &
        (Set([left_kink.k, left_kink.l]) == Set([right_kink.i, right_kink.j]))
    return true
  else
    return false
  end
end


"""Return pairs with indices of 'neighbouring' kinks that are Type-B-Entangled.
This function has no use in the current update set. Removability of a Kink will therefore not be considered here.
'neighbouring' refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first. (for type-B-entanglement that always imples the
vice versa case)
The Set consists of the pairs where the Type-B-entanglement is oriented
to the left of the first τ."""
function left_type_B_pairs(ck)
    pairs_left = Set{Tuple{Int,Int}}()
    for (i,(τ,kink)) in enumerate(ck)
        i_left = prev_affecting(ck, orbs(kink), τ)
        if is_type_B(last(ck[τ_left]), kink)
            push!(pairs_left, (i, i_left))
        end
    end
    return pairs_left
end


"""Return pairs with indices of 'neighbouring' kinks that are Type-B-Entangeld. Removablility of Kinks will
no be looked at. 'neighbouring' refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first. (for type-B-entanglement that always imples the
vice versa case)
The Set consists of the pairs where the Type-B-entanglement is oriented
#to the right of the first τ."""
function right_type_B_pairs(ck)
    pairs_right = Set{Tuple{Int,Int}}()
    for (i,(τ,kink)) in enumerate(ck)
        i_right = next_affecting(ck, orbs(kink), τ)
        if is_type_B(kink, last(ck[i_right]))
            push!(pairs_right, (i, i_right))
        end
    end
    return pairs_right
end


"""Return pairs with indices of 'neighbouring' Kinks that are Type-B-Entangeld AND removable.
'neighbouring' refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the other and is removable stands first. (for type-B-entanglement that always imples the
vice versa case)
The Set consists of the pairs where the Type-B-entanglement is oriented
#to the right of the first τ."""
function right_type_B_removable_pairs(ck)
    pairs_right = Set{Tuple{Int,Int}}()
    for (i,(τ,kink)) in enumerate(ck)
        i_right = next_affecting(ck, orbs(kink), τ)
        if is_type_B(kink, last(ck[i_right]))
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
    possible_new_orb_a(occs, orb_c, orb_d)
This function will find possibilites for the choice of the first annihilator for an add_type_B-update when the two creators are already selected.
To check possibility this will in particular check conservation laws and the occupation of the resulting fourth orb.
"""
function possible_new_orb_a(occs, orb_c, orb_d)
    filter(orb_a -> !in(find_fourth_orb_for_kink(orb_a, orb_c, orb_d),occs) && (orb_a != find_fourth_orb_for_kink(orb_a, orb_c, orb_d)),
                    setdiff!(sphere_with_same_spin(orb_c, dk = ex_radius), occs))
end

function add_type_B(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64, Step}
    #sampling propability
    prop_prob = 1
    #get first τ
    τ1 = ImgTime(rand())
    #do not allow two kinks at the same time
    while haskey(c.kinks, τ1)
        τ1 = ImgTime(rand())
    end
    occs = occupations_at(c, τ1)
    orb_c = rand(occs)
    prop_prob *= 1.0/e.N
    orb_d = rand(drop(occs, orb_c))
    prop_prob *= 1.0/(e.N-1)
    opportunities_orb_a = possible_new_orb_a(occs, orb_c, orb_d)
    if isempty(opportunities_orb_a)
        return 1.0, Step()
    end
    orb_a = rand(opportunities_orb_a)
    orb_b = find_fourth_orb_for_kink(orb_a, orb_c, orb_d)
    @assert !((orb_a == orb_d) | (orb_b == orb_c))
    @assert (!in(orb_b,occs) & (orb_a != orb_b))
    #We will change the proposal probability after we get τ2


    #get τ2
    borders = τ_borders(c.kinks, Set([orb_a,orb_b,orb_c,orb_d]),τ1)
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
    #If there are no kinks effecting any of the new kinks orbitals then anyone
    #of the two τs can be first τ
    if borders[2] == 1
        if rand() < 0.5
            firstτ = τ1
            lastτ = τ2
        else
            firstτ = τ2
            lastτ = τ1
        end
        prop_prob *= 0.5
    #If there are kinks which τ is firstτ depends on which is the
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
    #blocked over all permutations (see function “abs_Woffdiag_element”),
    #therefore the updates where we end up with the same kinks but start building
    #the kink with a different Excitation will result in a different order off indices
    #and therefore considered a different Update (to compensate that we use a factor ¼
    #in the off diagonal Matrix element contribution) .

    #However we have to consider the different ways off getting to the same
    #update by choosing the imaginary times in a different order.

    #Therefore we modify the proposal_probability in the following way
    occs_τ2 = occupations_at(c, τ2)
    opportunities_orb_a_τ2 = possible_new_orb_a(occs_τ2, orb_c, orb_d)
    @assert !isempty(opportunities_orb_a_τ2)
    prop_prob *= (1.0/length(opportunities_orb_a) + 1.0/length(opportunities_orb_a_τ2)) * 1.0/float(possible_τ2_interval)

    add_kinks = (
                firstτ => T4(orb_a,orb_b,orb_c,orb_d),
                # shuffle creators and shuffle annihilators
                lastτ => T4( random_shuffle(orb_d,orb_c)..., random_shuffle(orb_b, orb_a)... )
                )
    prop_prob *= 1/4

    #Look at wether the new pair of kinks modifies the Occupations at τ = 0
    if firstτ > lastτ
        drop_orbs = Set([orb_c, orb_d])
        add_orbs = Set([orb_a, orb_b])
    else
        drop_orbs = nothing
        add_orbs = Set{basis(c)}()
    end

    # MC step generated by this update
    Δ = Step(drop_orbs, Configuration(add_orbs, add_kinks...))

    # quotient of proposal probabilities
    dv = ( 1.0/ length( right_type_B_removable_pairs(apply_step(c,Δ).kinks) ) ) / prop_prob

    #calculate change in diagonal interaction energy
    delta_di = ΔWdiag_element(m, e, c, orb_a, orb_b, orb_c, orb_d, firstτ, lastτ)
    # weight factor
    dw = ((e.β)^2) * abs(Woffdiag_element(m, e, orb_a,orb_b,orb_c,orb_d))^2 * exp(-((delta_τ)*e.β * (energy(m, orb_a) + energy(m, orb_b) - energy(m, orb_c) - energy(m, orb_d)) + delta_di*e.β))


    @assert (dv*dw) >= 0
    @assert !isinf(dv*dw)
    @assert !isnan(dv*dw)
    occupations_at(c, ImgTime(0.99))
    return dv*dw, Δ
end

function remove_type_B(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64, Step}
    if isempty(c.kinks)
        return 1.0, Step()
    end
    opportunities = right_type_B_removable_pairs(c.kinks)
    if isempty(opportunities)
        return 1.0, Step()
    end
    #If a kink1 is entangeld to the right with the nearest kink, that
    # acts on one of his orbs, in a type-B-way that implies the vice versa case
    # therefor there is no value in distingusihing between left and right entanglement
    i1, i2 = rand(opportunities)
    kink1 = c.kinks[i1]
    kink2 = c.kinks[i2]
    prop_prob = 1.0/length(opportunities)

    #change configuration
    #see if occupations at τ=0 are modified
    if first(kink1) > first(kink2)
        # change_occupations_at(c.occupations, last(kink2))
        add_orbs = Set([last(kink2).i, last(kink2).j])
        drop_orbs = Set([last(kink2).k, last(kink2).l])
        delta_τ = first(kink2)-first(kink1) + 1.0
    else
        add_orbs = nothing
        drop_orbs = Set{basis(c)}()
        delta_τ = first(kink2)-first(kink1)
    end

    # MC step generated by this update
    Δ = Step(Configuration(drop_orbs, kink1, kink2), add_orbs)
    @assert(delta_τ > 0)
    @assert(delta_τ <= 1)

    c_new = apply_step(c,Δ)

    # calculate inverse prop_prob (see add_type_B)
    ijkl = Set([last(kink1).i, last(kink1).j, last(kink1).k, last(kink1).l])
    borders = τ_borders(c_new.kinks, ijkl, first(kink1))
    possible_τ2_interval = borders[2]-borders[1]
    if possible_τ2_interval < 0
        possible_τ2_interval = 1 + possible_τ2_interval
    end
    occs_τ_kink1 = occupations_at(c_new, first(kink1))
    occs_τ_kink2 = occupations_at(c_new, first(kink2))
    orb_a = last(kink1).i
    orb_b = last(kink1).j
    orb_c = last(kink1).k
    orb_d = last(kink1).l

    #See how prop_prob changes in the function add_type_B to understand this expression
    inverse_prop_prob = (1/e.N)*(1/(e.N-1)) * (1.0/length(possible_new_orb_a(occs_τ_kink1, orb_c, orb_d)) + 1/length(possible_new_orb_a(occs_τ_kink2, orb_c, orb_d))) * 1.0/float(possible_τ2_interval) * (1/4)


    if borders[2] == 1
        inverse_prop_prob *= 0.5
    end
    # quotient of proposal probabilities
    dv = inverse_prop_prob/prop_prob

    # calculate change in diagonal interaction energy
    delta_di = ΔWdiag_element(m, e, c_new, last(kink1).i, last(kink1).j, last(kink1).k, last(kink1).l, first(kink1), first(kink2))
    # weight factor
    dw = (1.0/(e.β)^2) * (1.0/(abs(Woffdiag_element(m, e, last(kink1))))^2) * exp((delta_τ)*e.β * (energy(m, orb_a) + energy(m, orb_b) - energy(m, orb_c) - energy(m, orb_d)) + e.β * delta_di)


    @assert (dv*dw) >= 0
    @assert !isinf(dv*dw)
    return dv*dw, Δ
end


function change_type_B(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64,Step} #This update is redundant wif we have add- and remove-Type-C-Updates
    if isempty(c.kinks)
        return 1.0, Step()
    end

    kink_opportunities = right_type_B_pairs(c.kinks)
    if isempty(kink_opportunities)
        return 1.0, Step()
    end
    i1, i2 = rand(kink_opportunities)
    kink1 = c.kinks[i1]
    kink2 = c.kinks[i2]
    occs = occupations_at(c, first(kink1))


    opportunities = filter( x -> isunaffected_in_interval(c.kinks, x, first(kink1), first(kink2))
                        ,possible_new_orb_a(occs, last(kink1).i, last(kink1).j))
    delete!(opportunities, last(kink1).k)
    delete!(opportunities, last(kink1).l)

    if isempty(opportunities)
        return 1.0, Step()
    end

    new_orb_i = rand(opportunities)
    new_orb_j = find_fourth_orb_for_kink(new_orb_i, last(kink1).i, last(kink1).j)
    if new_orb_i == new_orb_j
        return 1.0, Step()
    end
    if !isunaffected_in_interval(c.kinks,new_orb_j,first(kink1),first(kink2)) | in(new_orb_j,occs)
        return 1.0, Step()
    end

    # see if occupations change
    if first(kink1) > first(kink2)
        drop_orbs = creators(last(kink1))
        add_orbs = (new_orb_i, new_orb_j)
        delta_τ = first(kink2)-first(kink1) + 1.0
    else
        add_orbs = ()
        drop_orbs = ()
        delta_τ = first(kink2)-first(kink1)
    end


    #calculate change in diagonal interaction energy
    delta_di = ΔWdiag_element(m, e, c, new_orb_i, new_orb_j, last(kink1).i, last(kink1).j, first(kink1), first(kink2))

    # shuffle indices of the new orbs in the second kink
    drop_kinks = (kink1,kink2)
    add_kinks = (
                first(kink1) => T4(new_orb_i, new_orb_j, last(kink1).k, last(kink1).l),
                # shuffle creators and shuffle annihilators
                first(kink2) => T4( random_shuffle(last(kink2).i, last(kink2).j)..., random_shuffle(new_orb_i, new_orb_j)...)
                )
    # do not add factor 1/4 to proposal porbability (prop_prob *= 1/4)
    # since this cancels with the reverse update in the quotient of proposal propabilities

    excite!(occs, new_orb_i, new_orb_j, last(kink1).i, last(kink1).j)

    # MC Step generated by this update
    Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))

    new_kinks = add(drop(c.kinks, drop_kinks), add_kinks)

    kink_opportunities_reverse = right_type_B_pairs(new_kinks)
    opportunities_reverse = filter( x -> isunaffected_in_interval(new_kinks, x, first(kink1), first(kink2))
                                ,possible_new_orb_a(occs, new_orb_i, new_orb_j))
    delete!(opportunities_reverse, last(kink1).k)
    delete!(opportunities_reverse, last(kink1).l)

    # quotient of proposal probabilities
    dv = length(opportunities)*length(kink_opportunities)/(length(opportunities_reverse)*length(kink_opportunities_reverse))

    # weight factor
    dw = exp(-(e.β*delta_τ*(energy(m, new_orb_i) + energy(m, new_orb_j) - energy(m, last(kink1).i) - energy(m, last(kink1).j)) + e.β*delta_di)) * ( Woffdiag_element(m, e, new_orb_i, new_orb_j, last(kink1).k, last(kink1).l ) / Woffdiag_element(m, e, last(kink1)) )^2

    @assert (dw * dv) >= 0
    @assert !isinf(dw * dv)
    return dw*dv, Δ
end

function isuseful(c::Configuration, up::typeof(remove_type_B))
    if isempty(c.kinks) || length(c.kinks) == 3
        return false
    else
        return true
    end
end

function isuseful(c::Configuration, up::typeof(change_type_B))
    if isempty(c.kinks) || length(c.kinks) == 3
        return false
    else
        return true
    end
end
