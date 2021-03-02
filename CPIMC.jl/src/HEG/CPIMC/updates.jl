const ex_radius = 3 #max Radius for exitation

function move_particle(c::Configuration, e::Ensemble)
    free_orbitals = get_non_interacting_orbs_of_set(c, c.occupations)
    if length(free_orbitals) == 0
        return 1
    else
        x = rand(get_non_interacting_orbs_of_set(c, c.occupations))
    end
    oe = get_non_interacting_orbs_of_set(c,setdiff!(get_sphere_with_same_spin(x, dk = ex_radius), c.occupations))

    #if there are no empty non interacting orbitals in neighbourhood make no change
    if length(oe) == 0
        return(1)
    end
    y = rand(oe)
    @assert x != y "same Configuration proposed."

    delta_di = get_change_diagonal_interaction(c, e, T2(y,x), ImgTime(0), ImgTime(1))
    @assert delta_di != Inf
    # weight change
    dw = exp(-(e.β*(get_energy(y)-get_energy(x)) + delta_di))

    # change Configuration
    delete!(c.occupations,x)
    push!(c.occupations,y)

    # get orbitals for reverse update
    oe2 = get_non_interacting_orbs_of_set(c,setdiff!(get_sphere_with_same_spin(y, dk = ex_radius), c.occupations))

    # quotient of proposal probabilities
    dv = length(oe)/length(oe2)
    @assert delta_di != Inf
    @assert (dv*dw) >= 0
    return (dv*dw)
end

function add_type_B(c::Configuration, e::Ensemble)
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
        return 1
    end
    orb_a = rand(opportunities_orb_a)
    orb_b = OrbitalHEG((orb_c.vec-orb_a.vec) + orb_d.vec,orb_d.spin)
    @assert !((orb_a == orb_d) | (orb_b == orb_c))
    if (!in(orb_b,opportunities_orb_b) | (orb_a == orb_b))
        return 1
    end

    #We will change the proposal probability after we get τ2


    #get τ2
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





    #calculate change in diagonal interaction energy
    delta_di = get_change_diagonal_interaction(c, e, T4(orb_a,orb_b,orb_c,orb_d), firstτ, lastτ)

    #change configuration
    c.kinks[firstτ] = T4(orb_a,orb_b,orb_c,orb_d)
    #shuffle index order of the second Kink
    if rand() < 0.5
        if rand() < 0.5
            #swap both
            c.kinks[lastτ] = T4(orb_d,orb_c,orb_b,orb_a)
        else
            #swap anihilators
            c.kinks[lastτ] = T4(orb_c,orb_d,orb_b,orb_a)
        end
    else
        if rand() < 0.5
            #swap creators
            c.kinks[lastτ] = T4(orb_d,orb_c,orb_a,orb_b)
        else
            #swap none
            c.kinks[lastτ] = T4(orb_c,orb_d,orb_a,orb_b)
        end
    end
    prop_prob *= 1/4

    #Look at wether the new pair of Kinks modifies the Occupations at τ = 0
    if firstτ > lastτ
        change_occupations(c.occupations, T4(orb_a,orb_b,orb_c,orb_d))
    end

    # quotient of proposal probabilities
    dv = (1/length(get_right_type_B_pairs(c)))/prop_prob

    # weight factor
    dw = ((e.β)^2) *
            get_abs_offdiagonal_element(e,c,T4(orb_a,orb_b,orb_c,orb_d))^2 *
            exp(-((delta_τ)*e.β * (get_energy(orb_a) +
                 get_energy(orb_b) - get_energy(orb_c) - get_energy(orb_d)) + delta_di))
    @assert (dv*dw) >= 0
    @assert (dv*dw) != Inf
    @assert !isnan(dv*dw)
    return(dv*dw)
end

function remove_type_B(c::Configuration, e::Ensemble)
    if length(c.kinks) == 0
        return 1
    end
    opportunities = get_right_type_B_pairs(c)
    if length(opportunities) == 0
        return(1)
    end
    #If a kink1 is entangeld to the right with the nearest kink, that
    # acts on one of his orbs, in a type-B-way that implies the vice versa case
    # therefor there is no value in distingusihing between left and right entanglement
    τ_kink1, τ_kink2 = rand(get_right_type_B_pairs(c))
    kink1 = (τ_kink1, c.kinks[τ_kink1])
    kink2 = (τ_kink2, c.kinks[τ_kink2])

    if last(kink1).i.spin != last(kink1).k.spin
        return 1
    end
    # if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it
    if dot(last(kink1).i.vec-last(kink1).k.vec, last(kink1).i.vec-last(kink1).k.vec) > ex_radius^2
        return 1
    end
    prop_prob = 1/length(opportunities)
    #change configuration
    delete!(c.kinks, first(kink1))
    delete!(c.kinks, first(kink2))
    #see if occupations at τ=0 are modified
    if first(kink1) > first(kink2)
        change_occupations(c.occupations, last(kink2))
        delta_τ = first(kink2)-first(kink1) + 1
    else
        delta_τ = first(kink2)-first(kink1)
    end
    @assert(delta_τ > 0)
    @assert(delta_τ <= 1)
    #calculate inverse prop_prob (see  add_type_B)
    ijkl = Set([last(kink1).i, last(kink1).j, last(kink1).k, last(kink1).l])
    borders = get_τ_borders(c, ijkl ,first(kink1))
    possible_τ2_interval = borders[2]-borders[1]
    if possible_τ2_interval < 0
        possible_τ2_interval = 1 + possible_τ2_interval
    end
    occs_τ_kink1 = get_occupations_at(c, first(kink1))
    occs_τ_kink2 = get_occupations_at(c, first(kink2))
    orb_a = last(kink1).i
    orb_b = last(kink1).j
    orb_c = last(kink1).k
    orb_d = last(kink1).l
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
    delta_di = get_change_diagonal_interaction(c, e, last(kink1), first(kink1), first(kink2))


    # weight factor
    dw = (1.0/(e.β)^2) *
        (1.0/(get_abs_offdiagonal_element(e,c,last(kink1)))^2) *
            exp((delta_τ)*e.β * (get_energy(orb_a) +
                 get_energy(orb_b) - get_energy(orb_c) - get_energy(orb_d)) + delta_di)
    @assert (dv*dw) >= 0
    @assert (dv*dw) != Inf
    return (dv*dw)
end


function change_type_B(c::Configuration, e::Ensemble) #This update is redundant wif we have add- and remove-Type-C-Updates
    if length(c.kinks) == 0
        return 1
    end

    kink_opportunities = get_right_type_B_pairs(c)
    if length(kink_opportunities) == 0
        return(1)
    end
    τ1, τ2 = rand(kink_opportunities)
    kink1 = (τ1, c.kinks[τ1])
    kink2 = (τ2, c.kinks[τ2])
    occs = get_occupations_at(c, first(kink1))

    opportunities = get_non_interacting_orbs_of_set_in_interval(
                        c,setdiff!(
                            get_sphere_with_same_spin(last(kink1).i, dk = ex_radius
                            ), occs
                        ),first(kink1),first(kink2)
                    )
    delete!(opportunities, last(kink1).k)
    delete!(opportunities, last(kink1).l)
    if length(opportunities) == 0
        return 1
    end
    new_orb_i = rand(opportunities)
    new_orb_j = OrbitalHEG(last(kink1).j.vec + last(kink1).i.vec - new_orb_i.vec, last(kink1).j.spin)
    if new_orb_i == new_orb_j
        return 1
    end
    if (!is_non_interacting_in_interval(c,new_orb_j,first(kink1),first(kink2)) |
        in(new_orb_j,occs))
        return 1
    else
        #calculate change in diagonal interaction energy
        delta_di = get_change_diagonal_interaction(c, e, T4(new_orb_i, new_orb_j, last(kink1).i, last(kink1).j), first(kink1), first(kink2))

        #change occupations
        if first(kink1) > first(kink2)
            change_occupations(c.occupations, T4(new_orb_i, new_orb_j, last(kink1).i, last(kink1).j))
            delta_τ = first(kink2)-first(kink1) + 1
        else
            delta_τ = first(kink2)-first(kink1)
        end


        dw = exp(-(e.β*delta_τ*(get_energy(new_orb_i) + get_energy(new_orb_j) -
                                    get_energy(last(kink1).i) - get_energy(last(kink1).j)) + delta_di)) *
            (get_abs_offdiagonal_element(e,c,T4(new_orb_i, new_orb_j, last(kink1).k, last(kink1).l))/
                        get_abs_offdiagonal_element(e,c,last(kink1)))^2

        #change Kinks

        c.kinks[first(kink1)] = T4(new_orb_i, new_orb_j, last(kink1).k, last(kink1).l)
        if rand() <= 0.5
            c.kinks[first(kink2)] = T4(last(kink2).i, last(kink2).j, new_orb_i, new_orb_j,)
        else
            c.kinks[first(kink2)] = T4(last(kink2).i, last(kink2).j, new_orb_j, new_orb_i,)
        end
        change_occupations(occs, T4(new_orb_i, new_orb_j, last(kink1).i, last(kink1).j))
        opportunites_reverse = get_non_interacting_orbs_of_set_in_interval(
                                    c,setdiff!(
                                        get_sphere_with_same_spin(new_orb_i, dk = ex_radius
                                        ), occs
                                    ),first(kink1),first(kink2)
                                )
        delete!(opportunites_reverse, last(kink1).k)
        delete!(opportunites_reverse, last(kink1).l)
        @assert (dw * length(opportunities)/length(opportunites_reverse)) >= 0
        return (dw * length(opportunities)/length(opportunites_reverse))
    end
end


function add_type_C(c::Configuration, e::Ensemble)
    prop_prob = 1
    if length(c.kinks) == 0
        return 1
    end
    old_kink = rand(c.kinks)
    prop_prob *= 1/length(c.kinks)
    occs = get_occupations_at(c, first(old_kink))
    prop_prob *= 0.5 #left or right
    if rand() >= 0.5
        #add kink left
        opportunities_new_orb1 = setdiff!(get_sphere_with_same_spin(last(old_kink).k, dk = ex_radius), occs)
        delete!(opportunities_new_orb1, last(old_kink).k)
        delete!(opportunities_new_orb1, last(old_kink).l)
        if length(opportunities_new_orb1) == 0
            return 1
        end
        new_orb1 = rand(opportunities_new_orb1)
        prop_prob *= 1/length(opportunities_new_orb1)
        new_orb2 = OrbitalHEG(last(old_kink).j.vec + (last(old_kink).i.vec - new_orb1.vec), last(old_kink).l.spin)
        if in(new_orb2, occs) | (new_orb1 == new_orb2)
            return 1
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

        prop_prob *= 1/Float64(τ_Intervall)

        delta_di = get_change_diagonal_interaction(c, e, T4(new_orb1,new_orb2,last(old_kink).k,last(old_kink).l), τ_new_kink, first(old_kink))

        #change_Configuration
        #see if c.occupations change
        if τ_new_kink > first(old_kink)  #new kink was added left of old kink
            change_occupations(c.occupations, T4(new_orb1,new_orb2,last(old_kink).k,last(old_kink).l))
        end
        c.kinks[τ_new_kink] = T4(new_orb1, new_orb2, last(old_kink).k, last(old_kink).l)
        c.kinks[first(old_kink)] = T4(last(old_kink).i, last(old_kink).j, new_orb1, new_orb2)

        #calculate weight differance
        dw_off_diag = get_abs_offdiagonal_element(e,c,c.kinks[τ_new_kink]) * get_abs_offdiagonal_element(e,c,c.kinks[first(old_kink)]) /
                                                get_abs_offdiagonal_element(e,c,last(old_kink))
        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(get_energy(new_orb1) + get_energy(new_orb2) -
                                    get_energy(last(old_kink).k) - get_energy(last(old_kink).l)) + delta_di))

        inverse_prop_prob = (1/length(get_right_type_C_pairs(c))) * 0.5

        #shuffle annihilators of the changed kink
        if rand() < 0.5
            c.kinks[first(first(old_kink))] = T4(c.kinks[first(first(old_kink))].i, c.kinks[first(first(old_kink))].j,
                                                       c.kinks[first(first(old_kink))].l, c.kinks[first(first(old_kink))].k)
        else
            nothing
        end
        prop_prob *= 0.5
    else
        #add kink right
        opportunities_new_orb1 = setdiff!(get_sphere_with_same_spin(last(old_kink).i, dk = ex_radius), occs)
        delete!(opportunities_new_orb1, last(old_kink).k)
        delete!(opportunities_new_orb1, last(old_kink).l)
        if length(opportunities_new_orb1) == 0
            return 1
        end
        new_orb1 = rand(opportunities_new_orb1)
        prop_prob *= 1/length(opportunities_new_orb1)
        new_orb2 = OrbitalHEG(last(old_kink).l.vec + (last(old_kink).k.vec - new_orb1.vec), last(old_kink).j.spin)
        if in(new_orb2, occs) | (new_orb1 == new_orb2)
            return 1
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

        prop_prob *= 1/Float64(τ_Intervall)

        delta_di = get_change_diagonal_interaction(c, e, T4(new_orb1,new_orb2,last(old_kink).i,last(old_kink).j), first(old_kink), τ_new_kink)

        #change_Configuration
        #see if c.occupations change
        if τ_new_kink < first(old_kink)  #new kink was added right of old kink
            change_occupations(c.occupations, T4(new_orb1,new_orb2,last(old_kink).i,last(old_kink).j))
        end
        c.kinks[τ_new_kink] = T4(last(old_kink).i, last(old_kink).j, new_orb1, new_orb2)
        c.kinks[first(old_kink)] = T4(new_orb1, new_orb2, last(old_kink).k, last(old_kink).l)



        dw_off_diag = get_abs_offdiagonal_element(e,c,c.kinks[τ_new_kink]) * get_abs_offdiagonal_element(e,c,c.kinks[first(old_kink)]) /
                                                get_abs_offdiagonal_element(e,c,last(old_kink))

        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(get_energy(new_orb1) + get_energy(new_orb2) -
                                    get_energy(last(old_kink).i) - get_energy(last(old_kink).j)) + delta_di))

        inverse_prop_prob = (1/length(get_left_type_C_pairs(c))) * 0.5


        #shuffle creators of the changed kink
        if rand() < 0.5
            c.kinks[first(first(old_kink))] = T4(c.kinks[first(first(old_kink))].j, c.kinks[first(first(old_kink))].i,
                                                        c.kinks[first(first(old_kink))].k, c.kinks[first(first(old_kink))].l)
        else
            nothing
        end
        prop_prob *= 0.5
    end


    #check if sign was changend
    #to calculate the sign change use only orbitals old kink and the new orbs,
    #so swaping of ij or kl in one of the kinks does not effect the signchange
    signum = 1
    if dot((last(old_kink).j.vec - new_orb1.vec),(last(old_kink).j.vec - new_orb1.vec)) >
            dot((last(old_kink).i.vec - new_orb1.vec),(last(old_kink).i.vec - new_orb1.vec))
        signum*= -1
    end
    if dot((last(old_kink).l.vec - new_orb1.vec),(last(old_kink).l.vec - new_orb1.vec)) >
            dot((last(old_kink).k.vec - new_orb1.vec),(last(old_kink).k.vec - new_orb1.vec))
        signum*= -1
    end
    if dot((last(old_kink).i.vec - last(old_kink).l.vec),(last(old_kink).i.vec - last(old_kink).l.vec)) >
            dot((last(old_kink).i.vec - last(old_kink).k.vec),(last(old_kink).i.vec - last(old_kink).k.vec))
        signum*= -1
    end

    c.sign *= signum


    @assert(delta_τ > 0 )
    @assert(!isnan((inverse_prop_prob/prop_prob)*dw))
    return((inverse_prop_prob/prop_prob)*dw)
end

function remove_type_C(c::Configuration, e::Ensemble)
    prop_prob = 0.5
    if rand() > 0.5
        #removed kink left of changed kink
        opportunities = get_right_type_C_pairs(c)
        if length(opportunities) == 0
            return 1
        end
        removed_kink_τ, changed_kink_τ = rand(opportunities)
        prop_prob *= 1/length(opportunities)

        if c.kinks[removed_kink_τ].i.spin != c.kinks[removed_kink_τ].k.spin
            @assert false
            return 1
        end
        # if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it
        if dot(c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec,
                    c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec) > (ex_radius^2)
            return 1
        end


        #safe thoose for later
        removed_orb1 = c.kinks[changed_kink_τ].k
        removed_orb2 = c.kinks[changed_kink_τ].l

        #change configuration
        c.kinks[changed_kink_τ] = T4(c.kinks[changed_kink_τ].i, c.kinks[changed_kink_τ].j, c.kinks[removed_kink_τ].k,c.kinks[removed_kink_τ].l)
        @assert removed_orb1 != c.kinks[changed_kink_τ].k
        #see if c.occupations change
        if removed_kink_τ > changed_kink_τ
            change_occupations(c.occupations, T4(c.kinks[changed_kink_τ].k,c.kinks[changed_kink_τ].l,removed_orb1,removed_orb2))
        end

        delete!(c.kinks, removed_kink_τ)

        #calculate reverse_prop_prob
        occs = get_occupations_at(c, changed_kink_τ)
        opportunities_reverse_new_orb1 = setdiff!(get_sphere_with_same_spin(c.kinks[changed_kink_τ].k, dk = ex_radius), occs)
        delete!(opportunities_reverse_new_orb1, c.kinks[changed_kink_τ].k)
        delete!(opportunities_reverse_new_orb1, c.kinks[changed_kink_τ].l)

        τ_Intervall = changed_kink_τ - first(get_τ_borders(c, Set([removed_orb1, removed_orb2,
                         c.kinks[changed_kink_τ].k, c.kinks[changed_kink_τ].l]),changed_kink_τ))


        if τ_Intervall < 0
            τ_Intervall +=1
        end

        inverse_prop_prob = (0.5/length(c.kinks)) * (1/length(opportunities_reverse_new_orb1)) *
                                 (1/Float64(τ_Intervall)) * (1/2)

        #calculate weight change
        delta_di = get_change_diagonal_interaction(c, e, T4(removed_orb1,removed_orb2, c.kinks[changed_kink_τ].k,c.kinks[changed_kink_τ].l), removed_kink_τ, changed_kink_τ)

        dw_off_diag = get_abs_offdiagonal_element(e,c,T4(removed_orb1, removed_orb2,  c.kinks[changed_kink_τ].k, c.kinks[changed_kink_τ].l)) *
                        get_abs_offdiagonal_element(e,c,T4(c.kinks[changed_kink_τ].i, c.kinks[changed_kink_τ].j, removed_orb1, removed_orb2)) /
                            get_abs_offdiagonal_element(e,c,c.kinks[changed_kink_τ])

        delta_τ = Float64(changed_kink_τ - removed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1/e.β)* (1/dw_off_diag) * exp(e.β * delta_τ * (get_energy(removed_orb1) + get_energy(removed_orb2) -
                                    get_energy(c.kinks[changed_kink_τ].k) - get_energy(c.kinks[changed_kink_τ].l)) + delta_di)

    else
        #removed kink right of changed kink
        opportunities = get_left_type_C_pairs(c)
        if length(opportunities) == 0
            return 1
        end
        removed_kink_τ, changed_kink_τ = rand(opportunities)
        prop_prob *= 1/length(opportunities)

        if c.kinks[removed_kink_τ].i.spin != c.kinks[removed_kink_τ].k.spin
            return 1
        end
        # if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it
        if dot(c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec,
                    c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec) > (ex_radius^2)
            return 1
        end

        #safe thoose for later
        removed_orb1 = c.kinks[changed_kink_τ].i
        removed_orb2 = c.kinks[changed_kink_τ].j
        #change configuration
        c.kinks[changed_kink_τ] = T4(c.kinks[removed_kink_τ].i, c.kinks[removed_kink_τ].j, c.kinks[changed_kink_τ].k,c.kinks[changed_kink_τ].l)

        #see if c.occupations change
        if removed_kink_τ < changed_kink_τ  #new kink was added right of old kink
            change_occupations(c.occupations, T4(c.kinks[changed_kink_τ].i,c.kinks[changed_kink_τ].j,removed_orb1,removed_orb2))
        end

        delete!(c.kinks, removed_kink_τ)

        #calculate reverse_prop_prob
        occs = get_occupations_at(c, changed_kink_τ)
        opportunities_reverse_new_orb1 = setdiff!(get_sphere_with_same_spin(c.kinks[changed_kink_τ].i, dk = ex_radius), occs)
        delete!(opportunities_reverse_new_orb1, c.kinks[changed_kink_τ].k)
        delete!(opportunities_reverse_new_orb1, c.kinks[changed_kink_τ].l)

        τ_Intervall =  last(get_τ_borders(c, Set([c.kinks[changed_kink_τ].i, c.kinks[changed_kink_τ].j,
                                                    removed_orb1, removed_orb2]),changed_kink_τ)) - changed_kink_τ

        if τ_Intervall < 0
            τ_Intervall +=1
        end
        inverse_prop_prob = (0.5/length(c.kinks)) * (1/length(opportunities_reverse_new_orb1)) *
                                 (1/Float64(τ_Intervall)) * (1/2)

        #calculate weight change
        delta_di = get_change_diagonal_interaction(c, e, T4(removed_orb1,removed_orb2, c.kinks[changed_kink_τ].i,c.kinks[changed_kink_τ].j), changed_kink_τ, removed_kink_τ)

        dw_off_diag = get_abs_offdiagonal_element(e,c,T4(removed_orb1, removed_orb2,  c.kinks[changed_kink_τ].k, c.kinks[changed_kink_τ].l)) *
                        get_abs_offdiagonal_element(e,c,T4(c.kinks[changed_kink_τ].i, c.kinks[changed_kink_τ].j, removed_orb1, removed_orb2)) /
                            get_abs_offdiagonal_element(e,c,c.kinks[changed_kink_τ])

        delta_τ = Float64(removed_kink_τ - changed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1/e.β)*(1/dw_off_diag) * exp(e.β * delta_τ*(get_energy(removed_orb1) + get_energy(removed_orb2) -
                                    get_energy(c.kinks[changed_kink_τ].i) - get_energy(c.kinks[changed_kink_τ].j)) + delta_di)

    end
    #check if sign was changend
    signum = 1
    if dot((c.kinks[changed_kink_τ].j.vec - removed_orb1.vec),(c.kinks[changed_kink_τ].j.vec - removed_orb1.vec)) >
            dot((c.kinks[changed_kink_τ].i.vec - removed_orb1.vec),(c.kinks[changed_kink_τ].i.vec - removed_orb1.vec))
        signum*= -1
    end
    if dot((c.kinks[changed_kink_τ].l.vec - removed_orb1.vec),(c.kinks[changed_kink_τ].l.vec - removed_orb1.vec)) >
            dot((c.kinks[changed_kink_τ].k.vec - removed_orb1.vec),(c.kinks[changed_kink_τ].k.vec - removed_orb1.vec))
        signum*= -1
    end
    if dot((c.kinks[changed_kink_τ].i.vec -c.kinks[changed_kink_τ].l.vec),(c.kinks[changed_kink_τ].i.vec - c.kinks[changed_kink_τ].l.vec)) >
            dot((c.kinks[changed_kink_τ].i.vec - c.kinks[changed_kink_τ].k.vec),(c.kinks[changed_kink_τ].i.vec - c.kinks[changed_kink_τ].k.vec))
        signum*= -1
    end
    c.sign *= signum


    @assert(delta_τ > 0 )
    @assert(!isnan((inverse_prop_prob/prop_prob) * dw))
    return((inverse_prop_prob/prop_prob) * dw)
end





function add_type_D(c::Configuration, e::Ensemble)
    prop_prob = 1
    if length(c.kinks) == 0
        return 1
    end
    old_kink = rand(c.kinks)
    prop_prob *= 1/length(c.kinks)
    occs = get_occupations_at(c, first(old_kink))
    prop_prob *= 0.5 #left or right
    if rand() >= 0.5
        #add kink left
        opportunities_new_orb1 = intersect!(get_sphere_with_same_spin(last(old_kink).i, dk = ex_radius), occs)
        delete!(opportunities_new_orb1, last(old_kink).i)
        delete!(opportunities_new_orb1, last(old_kink).j)
        if length(opportunities_new_orb1) == 0
            return 1
        end
        new_orb1 = rand(opportunities_new_orb1)
        prop_prob *= 1/length(opportunities_new_orb1)
        new_orb2 = OrbitalHEG(last(old_kink).k.vec + (last(old_kink).l.vec - new_orb1.vec), last(old_kink).j.spin)
        if !in(new_orb2, occs) | (new_orb1 == new_orb2)
            return 1
        end
        τ_Intervall = first(old_kink) - first(get_τ_borders(c, Set([last(old_kink).i, last(old_kink).j, new_orb1, new_orb2]), first(old_kink)))

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

        prop_prob *= 1/Float64(τ_Intervall)

        delta_di = get_change_diagonal_interaction(c, e, T4(last(old_kink).i,last(old_kink).j, new_orb1, new_orb2), τ_new_kink, first(old_kink))

        #change_Configuration
        #see if c.occupations change
        if τ_new_kink > first(old_kink)  #new kink was added left of old kink
            change_occupations(c.occupations, T4(last(old_kink).i,last(old_kink).j, new_orb1, new_orb2))
        end
        c.kinks[τ_new_kink] = T4(last(old_kink).i, last(old_kink).j,new_orb1, new_orb2)
        c.kinks[first(old_kink)] = T4(new_orb1, new_orb2, last(old_kink).k, last(old_kink).l)

        #calculate weight differance
        dw_off_diag = get_abs_offdiagonal_element(e,c,c.kinks[τ_new_kink]) * get_abs_offdiagonal_element(e,c,c.kinks[first(old_kink)]) /
                                                get_abs_offdiagonal_element(e,c,last(old_kink))
        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(get_energy(last(old_kink).i) + get_energy(last(old_kink).j)-
                                                         get_energy(new_orb1) - get_energy(new_orb2)) + delta_di))

        inverse_prop_prob = (1/length(get_right_type_D_pairs(c))) * 0.5

        #shuffle annihilators of the changed kink
        if rand() < 0.5
            c.kinks[first(first(old_kink))] = T4(c.kinks[first(first(old_kink))].j, c.kinks[first(first(old_kink))].i,
                                                       c.kinks[first(first(old_kink))].k, c.kinks[first(first(old_kink))].l)
        else
            nothing
        end
        prop_prob *= 0.5
    else
        #add kink right
        opportunities_new_orb1 = intersect!(get_sphere_with_same_spin(last(old_kink).k, dk = ex_radius), occs)
        delete!(opportunities_new_orb1, last(old_kink).i)
        delete!(opportunities_new_orb1, last(old_kink).j)
        if length(opportunities_new_orb1) == 0
            return 1
        end
        new_orb1 = rand(opportunities_new_orb1)
        prop_prob *= 1/length(opportunities_new_orb1)
        new_orb2 = OrbitalHEG(last(old_kink).i.vec + (last(old_kink).j.vec - new_orb1.vec), last(old_kink).l.spin)
        if !in(new_orb2, occs) | (new_orb1 == new_orb2)
            return 1
        end
        τ_Intervall = last(get_τ_borders(c, Set([
                        last(old_kink).k, last(old_kink).l, new_orb1, new_orb2]),first(old_kink))) -
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

        prop_prob *= 1/Float64(τ_Intervall)

        delta_di = get_change_diagonal_interaction(c, e, T4(last(old_kink).k, last(old_kink).l,new_orb1,new_orb2), first(old_kink), τ_new_kink)

        #change_Configuration
        #see if c.occupations change
        if τ_new_kink < first(old_kink)  #new kink was added right of old kink
            change_occupations(c.occupations, T4(last(old_kink).k, last(old_kink).l,new_orb1,new_orb2))
        end
        c.kinks[τ_new_kink] = T4(new_orb1, new_orb2, last(old_kink).k, last(old_kink).l)
        c.kinks[first(old_kink)] = T4(last(old_kink).i, last(old_kink).j, new_orb1, new_orb2)



        dw_off_diag = get_abs_offdiagonal_element(e,c,c.kinks[τ_new_kink]) * get_abs_offdiagonal_element(e,c,c.kinks[first(old_kink)]) /
                                                get_abs_offdiagonal_element(e,c,last(old_kink))

        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(get_energy(last(old_kink).k) + get_energy(last(old_kink).l) -
                                                         get_energy(new_orb1) - get_energy(new_orb2)) + delta_di))

        inverse_prop_prob = (1/length(get_left_type_D_pairs(c))) * 0.5


        #shuffle creators of the changed kink
        if rand() < 0.5
            c.kinks[first(first(old_kink))] = T4(c.kinks[first(first(old_kink))].i, c.kinks[first(first(old_kink))].j,
                                                        c.kinks[first(first(old_kink))].l, c.kinks[first(first(old_kink))].k)
        else
            nothing
        end
        prop_prob *= 0.5
    end


    #check if sign was changend
    signum = 1
    if dot((last(old_kink).j.vec - new_orb1.vec),(last(old_kink).j.vec - new_orb1.vec)) >
            dot((last(old_kink).i.vec - new_orb1.vec),(last(old_kink).i.vec - new_orb1.vec))
        signum*= -1
    end
    if dot((last(old_kink).l.vec - new_orb1.vec),(last(old_kink).l.vec - new_orb1.vec)) >
            dot((last(old_kink).k.vec - new_orb1.vec),(last(old_kink).k.vec - new_orb1.vec))
        signum*= -1
    end
    if dot((last(old_kink).i.vec - last(old_kink).l.vec),(last(old_kink).i.vec - last(old_kink).l.vec)) >
            dot((last(old_kink).i.vec - last(old_kink).k.vec),(last(old_kink).i.vec - last(old_kink).k.vec))
        signum*= -1
    end

    c.sign *= signum


    @assert(delta_τ > 0 )
    @assert(!isnan((inverse_prop_prob/prop_prob)*dw))
    return((inverse_prop_prob/prop_prob)*dw)
end

function remove_type_D(c::Configuration, e::Ensemble)
    prop_prob = 0.5
    if rand() > 0.5
        #removed kink left of changed kink
        opportunities = get_right_type_D_pairs(c)
        if length(opportunities) == 0
            return 1
        end
        removed_kink_τ, changed_kink_τ = rand(opportunities)
        prop_prob *= 1/length(opportunities)

        if c.kinks[removed_kink_τ].i.spin != c.kinks[removed_kink_τ].k.spin
            @assert false
            return 1
        end
        # if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it
        if dot(c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec,
                    c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec) > (ex_radius^2)
            return 1
        end


        #safe thoose for later
        removed_orb1 = c.kinks[removed_kink_τ].k
        removed_orb2 = c.kinks[removed_kink_τ].l

        #change configuration
        c.kinks[changed_kink_τ] = T4(c.kinks[removed_kink_τ].i, c.kinks[removed_kink_τ].j, c.kinks[changed_kink_τ].k,c.kinks[changed_kink_τ].l)
        @assert removed_orb1 != c.kinks[changed_kink_τ].k
        #see if c.occupations change
        if removed_kink_τ > changed_kink_τ
            change_occupations(c.occupations, T4(removed_orb1,removed_orb2, c.kinks[changed_kink_τ].i,c.kinks[changed_kink_τ].j))
        end

        delete!(c.kinks, removed_kink_τ)

        #calculate reverse_prop_prob
        occs = get_occupations_at(c, changed_kink_τ)
        opportunities_reverse_new_orb1 = intersect!(get_sphere_with_same_spin(c.kinks[changed_kink_τ].i, dk = ex_radius), occs)
        delete!(opportunities_reverse_new_orb1, c.kinks[changed_kink_τ].i)
        delete!(opportunities_reverse_new_orb1, c.kinks[changed_kink_τ].j)
        @assert(in(removed_orb1,opportunities_reverse_new_orb1))
        τ_Intervall = changed_kink_τ - first(get_τ_borders(c, Set([removed_orb1, removed_orb2,
                         c.kinks[changed_kink_τ].i, c.kinks[changed_kink_τ].j]),changed_kink_τ))


        if τ_Intervall < 0
            τ_Intervall +=1
        end

        inverse_prop_prob = (0.5/length(c.kinks)) * (1/length(opportunities_reverse_new_orb1)) *
                                 (1/Float64(τ_Intervall)) * (1/2)

        #calculate weight change
        delta_di = get_change_diagonal_interaction(c, e, T4(c.kinks[changed_kink_τ].i,c.kinks[changed_kink_τ].j, removed_orb1,removed_orb2), removed_kink_τ, changed_kink_τ)

        dw_off_diag = get_abs_offdiagonal_element(e,c,T4(removed_orb1, removed_orb2,  c.kinks[changed_kink_τ].k, c.kinks[changed_kink_τ].l)) *
                        get_abs_offdiagonal_element(e,c,T4(c.kinks[changed_kink_τ].i, c.kinks[changed_kink_τ].j, removed_orb1, removed_orb2)) /
                            get_abs_offdiagonal_element(e,c,c.kinks[changed_kink_τ])

        delta_τ = Float64(changed_kink_τ - removed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1/e.β)* (1/dw_off_diag) * exp(e.β * delta_τ * (get_energy(c.kinks[changed_kink_τ].i) + get_energy(c.kinks[changed_kink_τ].j) -
                                                                    get_energy(removed_orb1) - get_energy(removed_orb2)) + delta_di)

    else
        #removed kink right of changed kink
        opportunities = get_left_type_D_pairs(c)
        if length(opportunities) == 0
            return 1
        end
        removed_kink_τ, changed_kink_τ = rand(opportunities)
        prop_prob *= 1/length(opportunities)

        if c.kinks[removed_kink_τ].i.spin != c.kinks[removed_kink_τ].k.spin
            return 1
        end
        # if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it
        if dot(c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec,
                    c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec) > (ex_radius^2)
            return 1
        end

        #safe thoose for later
        removed_orb1 = c.kinks[removed_kink_τ].i
        removed_orb2 = c.kinks[removed_kink_τ].j
        #change configuration
        c.kinks[changed_kink_τ] = T4(c.kinks[changed_kink_τ].i, c.kinks[changed_kink_τ].j, c.kinks[removed_kink_τ].k,c.kinks[removed_kink_τ].l)

        #see if c.occupations change
        if removed_kink_τ < changed_kink_τ  #new kink was added right of old kink
            change_occupations(c.occupations, T4(removed_orb1, removed_orb2, c.kinks[changed_kink_τ].k,c.kinks[changed_kink_τ].l))
        end

        delete!(c.kinks, removed_kink_τ)

        #calculate reverse_prop_prob
        occs = get_occupations_at(c, changed_kink_τ)
        opportunities_reverse_new_orb1 = intersect!(get_sphere_with_same_spin(c.kinks[changed_kink_τ].k, dk = ex_radius), occs)
        delete!(opportunities_reverse_new_orb1, c.kinks[changed_kink_τ].i)
        delete!(opportunities_reverse_new_orb1, c.kinks[changed_kink_τ].j)
        @assert(in(removed_orb1,opportunities_reverse_new_orb1))
        τ_Intervall =  last(get_τ_borders(c, Set([c.kinks[changed_kink_τ].k, c.kinks[changed_kink_τ].l,
                                                    removed_orb1, removed_orb2]),changed_kink_τ)) - changed_kink_τ

        if τ_Intervall < 0
            τ_Intervall +=1
        end
        inverse_prop_prob = (0.5/length(c.kinks)) * (1/length(opportunities_reverse_new_orb1)) *
                                 (1/Float64(τ_Intervall)) * (1/2)

        #calculate weight change
        delta_di = get_change_diagonal_interaction(c, e, T4(c.kinks[changed_kink_τ].k,c.kinks[changed_kink_τ].l, removed_orb1,removed_orb2), changed_kink_τ, removed_kink_τ)

        dw_off_diag = get_abs_offdiagonal_element(e,c,T4(removed_orb1, removed_orb2,  c.kinks[changed_kink_τ].k, c.kinks[changed_kink_τ].l)) *
                        get_abs_offdiagonal_element(e,c,T4(c.kinks[changed_kink_τ].i, c.kinks[changed_kink_τ].j, removed_orb1, removed_orb2)) /
                            get_abs_offdiagonal_element(e,c,c.kinks[changed_kink_τ])

        delta_τ = Float64(removed_kink_τ - changed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1/e.β)*(1/dw_off_diag) * exp(e.β * delta_τ*(get_energy(c.kinks[changed_kink_τ].k) + get_energy(c.kinks[changed_kink_τ].l) -
                                                            get_energy(removed_orb1) - get_energy(removed_orb2)) + delta_di)

    end
    #check if sign was changend
    signum = 1
    if dot((c.kinks[changed_kink_τ].j.vec - removed_orb1.vec),(c.kinks[changed_kink_τ].j.vec - removed_orb1.vec)) >
            dot((c.kinks[changed_kink_τ].i.vec - removed_orb1.vec),(c.kinks[changed_kink_τ].i.vec - removed_orb1.vec))
        signum*= -1
    end
    if dot((c.kinks[changed_kink_τ].l.vec - removed_orb1.vec),(c.kinks[changed_kink_τ].l.vec - removed_orb1.vec)) >
            dot((c.kinks[changed_kink_τ].k.vec - removed_orb1.vec),(c.kinks[changed_kink_τ].k.vec - removed_orb1.vec))
        signum*= -1
    end
    if dot((c.kinks[changed_kink_τ].i.vec -c.kinks[changed_kink_τ].l.vec),(c.kinks[changed_kink_τ].i.vec - c.kinks[changed_kink_τ].l.vec)) >
            dot((c.kinks[changed_kink_τ].i.vec - c.kinks[changed_kink_τ].k.vec),(c.kinks[changed_kink_τ].i.vec - c.kinks[changed_kink_τ].k.vec))
        signum*= -1
    end
    c.sign *= signum

    @assert(dw != Inf)
    @assert(delta_τ > 0 )
    @assert(!isnan((inverse_prop_prob/prop_prob) * dw))

    return((inverse_prop_prob/prop_prob) * dw)
end



function add_type_E(c::Configuration, e::Ensemble)
    #After the updates the i and l komponents of both kinks will contain the old kink, wile the j and k components contain the old orbitals
    prop_prob = 1
    if length(c.kinks) == 0
        return 1
    end
    old_kink = rand(c.kinks)
    prop_prob *= 1/length(c.kinks)
    occs = get_occupations_at(c, first(old_kink))
    prop_prob *= 0.5 #left or right
    if rand() >= 0.5
        #add kink left
        #now choose orbitals of the old Kinks that shell also be part of the new Kink
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
        opportunities_new_kink_new_annihilator = intersect!(get_sphere_with_same_spin(new_kink_old_creator, dk = ex_radius), occs)
        delete!(opportunities_new_kink_new_annihilator, last(old_kink).i)
        delete!(opportunities_new_kink_new_annihilator, last(old_kink).j)
        if length(opportunities_new_kink_new_annihilator) == 0
            return 1
        end
        new_kink_new_annihilator = rand(opportunities_new_kink_new_annihilator)
        prop_prob *= 1/length(opportunities_new_kink_new_annihilator)
        #calculate new creator
        new_kink_new_creator = OrbitalHEG(new_kink_old_annihilator.vec + (new_kink_new_annihilator.vec - new_kink_old_creator.vec), new_kink_old_annihilator.spin)
        if (in(new_kink_new_creator, occs) | (new_kink_new_creator == last(old_kink).k) | (new_kink_new_creator == last(old_kink).l))
            return 1
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

        prop_prob *= 1/Float64(τ_Intervall)

        delta_di = get_change_diagonal_interaction(c, e, T4(new_kink_old_creator,new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator), τ_new_kink, first(old_kink))

        #change_Configuration
        #see if c.occupations change
        if τ_new_kink > first(old_kink)  #consider that new kink was added left of old kink
            @assert (!in(new_kink_new_creator, c.occupations))
            change_occupations(c.occupations, T4(new_kink_old_creator,new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator))
        end
        c.kinks[τ_new_kink] = T4(new_kink_old_creator,new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator)
        c.kinks[first(old_kink)] = T4(changed_kink_old_creator, new_kink_new_annihilator, new_kink_new_creator, changed_kink_old_annihilator)

        #calculate weight differance
        dw_off_diag = get_abs_offdiagonal_element(e,c,c.kinks[τ_new_kink]) * get_abs_offdiagonal_element(e,c,c.kinks[first(old_kink)]) /
                                                get_abs_offdiagonal_element(e,c,last(old_kink))
        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(get_energy(c.kinks[τ_new_kink].i) + get_energy(c.kinks[τ_new_kink].j)-
                                                         get_energy(c.kinks[τ_new_kink].k) - get_energy(c.kinks[τ_new_kink].l)) + delta_di))

        inverse_prop_prob = (1/length(get_right_type_E_removable_pairs(c))) * 0.5 * 0.25

        @assert(inverse_prop_prob != Inf)
    else
        #add kink right
        #now choose orbitals of the old Kinks that shell also be part of the new Kink
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
        opportunities_new_kink_new_annihilator = setdiff!(get_sphere_with_same_spin(new_kink_old_creator, dk = ex_radius), occs)
        delete!(opportunities_new_kink_new_annihilator, last(old_kink).k)
        delete!(opportunities_new_kink_new_annihilator, last(old_kink).l)
        if length(opportunities_new_kink_new_annihilator) == 0
            return 1
        end
        new_kink_new_annihilator = rand(opportunities_new_kink_new_annihilator)
        prop_prob *= 1/length(opportunities_new_kink_new_annihilator)
        new_kink_new_creator = OrbitalHEG(new_kink_old_annihilator.vec + (new_kink_new_annihilator.vec - new_kink_old_creator.vec), new_kink_old_annihilator.spin)
        if (!in(new_kink_new_creator, occs) | (new_kink_new_creator == last(old_kink).i) | (new_kink_new_creator == last(old_kink).j))
            return 1
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

        prop_prob *= 1/Float64(τ_Intervall)
                                                          #Inverse new Kink
        delta_di = get_change_diagonal_interaction(c, e, T4(new_kink_new_annihilator, new_kink_old_annihilator, new_kink_old_creator, new_kink_new_creator), first(old_kink), τ_new_kink)

        #change_Configuration
        #see if c.occupations change
        if τ_new_kink < first(old_kink)  #new kink was added right of old kink
            change_occupations(c.occupations, T4(new_kink_new_annihilator, new_kink_old_annihilator, new_kink_old_creator, new_kink_new_creator))
        end
        c.kinks[τ_new_kink] = T4(new_kink_old_creator, new_kink_new_creator, new_kink_new_annihilator, new_kink_old_annihilator)
        c.kinks[first(old_kink)] = T4(changed_kink_old_creator, new_kink_new_annihilator, new_kink_new_creator, changed_kink_old_annihilator)



        dw_off_diag = get_abs_offdiagonal_element(e,c,c.kinks[τ_new_kink]) * get_abs_offdiagonal_element(e,c,c.kinks[first(old_kink)]) /
                                                get_abs_offdiagonal_element(e,c,last(old_kink))

        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(get_energy(c.kinks[τ_new_kink].k) + get_energy(c.kinks[τ_new_kink].l) -
                                                         get_energy(c.kinks[τ_new_kink].i) - get_energy(c.kinks[τ_new_kink].j)) + delta_di))

        inverse_prop_prob = (1/length(get_left_type_E_removable_pairs(c))) * 0.5 * 0.25


        @assert(inverse_prop_prob != Inf)
    end

    #check if sign was changend
    signum = 1
    if dot((new_kink_old_creator.vec - new_kink_new_annihilator.vec),(new_kink_old_creator.vec - new_kink_new_annihilator.vec)) >
            dot((new_kink_old_creator.vec - new_kink_old_annihilator.vec),(new_kink_old_creator.vec - new_kink_old_annihilator.vec))
        signum*= -1
    end
    if dot((changed_kink_old_creator.vec - new_kink_new_creator.vec),(changed_kink_old_creator.vec - new_kink_new_creator.vec)) >
            dot((changed_kink_old_creator.vec - changed_kink_old_annihilator.vec),(changed_kink_old_creator.vec - changed_kink_old_annihilator.vec))
        signum*= -1
    end
    if dot((new_kink_old_creator.vec - changed_kink_old_annihilator.vec),(new_kink_old_creator.vec - changed_kink_old_annihilator.vec)) >
            dot((new_kink_old_creator.vec - new_kink_old_annihilator.vec),(new_kink_old_creator.vec - new_kink_old_annihilator.vec))
        signum*= -1
    end

    c.sign *= signum


    #shuffle Indices
    #shuffle changed kink
    if rand() < 0.5
        if rand() < 0.5
            c.kinks[first(old_kink)] = T4(c.kinks[first(old_kink)].j, c.kinks[first(old_kink)].i,
                                                        c.kinks[first(old_kink)].l, c.kinks[first(old_kink)].k)
        else
            c.kinks[first(old_kink)] = T4(c.kinks[first(old_kink)].i, c.kinks[first(old_kink)].j,
                                                        c.kinks[first(old_kink)].l, c.kinks[first(old_kink)].k)
        end
    else
        if rand() < 0.5
            c.kinks[first(old_kink)] = T4(c.kinks[first(old_kink)].j, c.kinks[first(old_kink)].i,
                                                        c.kinks[first(old_kink)].k, c.kinks[first(old_kink)].l)
        else
            c.kinks[first(old_kink)] = T4(c.kinks[first(old_kink)].i, c.kinks[first(old_kink)].j,
                                                        c.kinks[first(old_kink)].k, c.kinks[first(old_kink)].l)
        end
    end
    #shuffle new Kink
    if rand() < 0.5
        if rand() < 0.5
            c.kinks[τ_new_kink] = T4(c.kinks[τ_new_kink].j, c.kinks[τ_new_kink].i,
                                                        c.kinks[τ_new_kink].l, c.kinks[τ_new_kink].k)
        else
            c.kinks[τ_new_kink] = T4(c.kinks[τ_new_kink].i, c.kinks[τ_new_kink].j,
                                                        c.kinks[τ_new_kink].l, c.kinks[τ_new_kink].k)
        end
    else
        if rand() < 0.5
            c.kinks[τ_new_kink] = T4(c.kinks[τ_new_kink].j, c.kinks[τ_new_kink].i,
                                                        c.kinks[τ_new_kink].k, c.kinks[τ_new_kink].l)
        else
            c.kinks[τ_new_kink] = T4(c.kinks[τ_new_kink].i, c.kinks[τ_new_kink].j,
                                                        c.kinks[τ_new_kink].k, c.kinks[τ_new_kink].l)
        end
    end
    @assert(is_type_E(c.kinks[τ_new_kink], c.kinks[first(old_kink)]) != false)
    prop_prob *= 1/16

    @assert(delta_τ > 0 )
    @assert(dw != Inf)
    @assert(prop_prob != 0)
    @assert(inverse_prop_prob != Inf)
    @assert(inverse_prop_prob != 0)
    @assert(!isnan((inverse_prop_prob/prop_prob)*dw))
    return((inverse_prop_prob/prop_prob)*dw)
end

function remove_type_E(c::Configuration, e::Ensemble)
    prop_prob = 0.5
    if rand() > 0.5
        #removed kink left of changed kink
        opportunities = get_right_type_E_removable_pairs(c)
        if length(opportunities) == 0
            return 1
        end
        removed_kink_tuple, changed_kink_tuple = rand(opportunities)
        removed_kink_τ =  first(removed_kink_tuple)
        changed_kink_τ =  first(changed_kink_tuple)
        prop_prob *= 1/length(opportunities)
        #safe thoose for later
        removed_Kink = last(removed_kink_tuple)
        changed_Kink_old = last(changed_kink_tuple)
        if removed_Kink.i.spin != removed_Kink.k.spin
            @assert false
            return 1
        end
        # if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it
        if dot(removed_Kink.i.vec-removed_Kink.k.vec,
                    removed_Kink.i.vec-removed_Kink.k.vec) > (ex_radius^2)
            @assert false #This should be prevented by using get_right_type_E_removable_pairs
            return 1
        end




        #change configuration
        c.kinks[changed_kink_τ] = T4(removed_Kink.i, changed_Kink_old.i, removed_Kink.l,changed_Kink_old.l)
        #see if c.occupations change
        if removed_kink_τ > changed_kink_τ
            change_occupations(c.occupations, T4(removed_Kink.k,removed_Kink.l, removed_Kink.i,removed_Kink.j))
        end
        delete!(c.kinks, removed_kink_τ)

        #calculate reverse_prop_prob
        occs = get_occupations_at(c, changed_kink_τ)
        opportunities_occ_orb_E = intersect!(get_sphere_with_same_spin(removed_Kink.i, dk = ex_radius), occs)
        delete!(opportunities_occ_orb_E, c.kinks[changed_kink_τ].i)
        delete!(opportunities_occ_orb_E, c.kinks[changed_kink_τ].j)
        @assert(in(removed_Kink.k,opportunities_occ_orb_E))
        τ_Intervall = changed_kink_τ - first(get_τ_borders(c, Set([removed_Kink.k, removed_Kink.l,
                         removed_Kink.i, removed_Kink.j]),changed_kink_τ))


        if τ_Intervall < 0
            τ_Intervall +=1
        end

        inverse_prop_prob = (0.5/length(c.kinks)) * (1/length(opportunities_occ_orb_E)) *
                                 (1/Float64(τ_Intervall)) * (1/4) * (1/16)

        #calculate weight change
        delta_di = get_change_diagonal_interaction(c, e, removed_Kink, removed_kink_τ, changed_kink_τ)

        dw_off_diag = get_abs_offdiagonal_element(e,c,removed_Kink) *
                        get_abs_offdiagonal_element(e,c,changed_Kink_old) /
                            get_abs_offdiagonal_element(e,c,c.kinks[changed_kink_τ])

        delta_τ = Float64(changed_kink_τ - removed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1/e.β)* (1/dw_off_diag) * exp(e.β * delta_τ * (get_energy(removed_Kink.i) + get_energy(removed_Kink.j) -
                                                                    get_energy(removed_Kink.k) - get_energy(removed_Kink.l)) + delta_di)

    else
        #removed kink right of changed kink
        opportunities = get_left_type_E_removable_pairs(c)
        if length(opportunities) == 0
            return 1
        end
        removed_kink_tuple, changed_kink_tuple = rand(opportunities)
        removed_kink_τ = first(removed_kink_tuple)
        changed_kink_τ = first(changed_kink_tuple)
        prop_prob *= 1/length(opportunities)
        removed_Kink = last(removed_kink_tuple)
        changed_Kink_old = last(changed_kink_tuple)
        if removed_Kink.i.spin != removed_Kink.k.spin
            return 1
        end
        # if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it
        if dot(removed_Kink.i.vec-removed_Kink.k.vec,
                    removed_Kink.i.vec-removed_Kink.k.vec) > (ex_radius^2)
            @assert false #This should not happen when using get_left_type_E_removable_pairs(c)
            return 1
        end


        #change configuration
        c.kinks[changed_kink_τ] = T4(removed_Kink.i, changed_Kink_old.i, removed_Kink.l, changed_Kink_old.l)

        #see if c.occupations change
        if removed_kink_τ < changed_kink_τ  #new kink was added right of old kink
            change_occupations(c.occupations, T4(removed_Kink.i, removed_Kink.j, removed_Kink.k,removed_Kink.l))
        end

        delete!(c.kinks, removed_kink_τ)

        #calculate reverse_prop_prob
        occs = get_occupations_at(c, changed_kink_τ)
        opportunities_unocc_orb_E = setdiff!(get_sphere_with_same_spin(removed_Kink.i, dk = ex_radius), occs)
        delete!(opportunities_unocc_orb_E, c.kinks[changed_kink_τ].k)
        delete!(opportunities_unocc_orb_E, c.kinks[changed_kink_τ].l)
        @assert(in(removed_Kink.k,opportunities_unocc_orb_E))
        τ_Intervall =  last(get_τ_borders(c, Set([removed_Kink.i, removed_Kink.j,
                                                    removed_Kink.k, removed_Kink.l]),changed_kink_τ)) - changed_kink_τ

        if τ_Intervall < 0
            τ_Intervall +=1
        end
        inverse_prop_prob = (0.5/length(c.kinks)) * (1/length(opportunities_unocc_orb_E)) *
                                 (1/Float64(τ_Intervall)) * (1/4) * (1/16)

        #calculate weight change
        delta_di = get_change_diagonal_interaction(c, e, T4(removed_Kink.k,removed_Kink.l, removed_Kink.i,removed_Kink.j), changed_kink_τ, removed_kink_τ)

        dw_off_diag = get_abs_offdiagonal_element(e,c,removed_Kink) *
                        get_abs_offdiagonal_element(e,c,changed_Kink_old) /
                            get_abs_offdiagonal_element(e,c,c.kinks[changed_kink_τ])

        delta_τ = Float64(removed_kink_τ - changed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1/e.β)*(1/dw_off_diag) * exp(e.β * delta_τ*(get_energy(removed_Kink.k) + get_energy(removed_Kink.l) -
                                                            get_energy(removed_Kink.i) - get_energy(removed_Kink.j)) + delta_di)

    end

    #check if sign was changend
    signum = 1
    if dot((changed_Kink_old.i.vec - changed_Kink_old.k.vec),(changed_Kink_old.i.vec - changed_Kink_old.k.vec)) >
            dot((changed_Kink_old.i.vec - changed_Kink_old.l.vec),(changed_Kink_old.i.vec - changed_Kink_old.l.vec))
        signum*= -1
    end
    if dot((removed_Kink.i.vec - removed_Kink.k.vec),(removed_Kink.i.vec - removed_Kink.k.vec)) >
            dot((removed_Kink.i.vec - removed_Kink.l.vec),(removed_Kink.i.vec - removed_Kink.l.vec))
        signum*= -1
    end
    if dot((changed_Kink_old.i.vec - removed_Kink.l.vec),(changed_Kink_old.i.vec - removed_Kink.l.vec)) >
            dot((changed_Kink_old.i.vec - changed_Kink_old.l.vec),(changed_Kink_old.i.vec - changed_Kink_old.l.vec))
        signum*= -1
    end
    c.sign *= signum

    #shuffle Indices
    if rand() < 0.5
        if rand() < 0.5
            c.kinks[changed_kink_τ] = T4(c.kinks[changed_kink_τ].j, c.kinks[changed_kink_τ].i,
                                                        c.kinks[changed_kink_τ].l, c.kinks[changed_kink_τ].k)
        else
            c.kinks[changed_kink_τ] = T4(c.kinks[changed_kink_τ].i, c.kinks[changed_kink_τ].j,
                                                        c.kinks[changed_kink_τ].l, c.kinks[changed_kink_τ].k)
        end
    else
        if rand() < 0.5
            c.kinks[changed_kink_τ] = T4(c.kinks[changed_kink_τ].j, c.kinks[changed_kink_τ].i,
                                                        c.kinks[changed_kink_τ].k, c.kinks[changed_kink_τ].l)
        else
            c.kinks[changed_kink_τ] = T4(c.kinks[changed_kink_τ].i, c.kinks[changed_kink_τ].j,
                                                        c.kinks[changed_kink_τ].k, c.kinks[changed_kink_τ].l)
        end
    end
    prop_prob *= 0.25


    @assert(dw != Inf)
    @assert(delta_τ > 0 )
    @assert(!isnan((inverse_prop_prob/prop_prob) * dw))

    return((inverse_prop_prob/prop_prob) * dw)
end




function shuffle_indices(c::Configuration, e::Ensemble)
    if length(c.kinks) == 0
        return(1)
    end
    kink = rand(c.kinks)
    if rand() > 0.5
        c.kinks[first(kink)] = T4(last(kink).j,last(kink).i,last(kink).k,last(kink).l)
    else
        c.kinks[first(kink)] = T4(last(kink).i,last(kink).j,last(kink).l,last(kink).k)
    end
    return(1)
end


function add_remove_kink_chain(c::Configuration,e::Ensemble)
    if length(c.kinks) == 0
        return(1)
    end
    add_single_kinks = [add_type_C, add_type_D, add_type_E]
    remove_single_kinks = [remove_type_C, remove_type_D, remove_type_E]
    number_of_kinks = rand([1,2,3,4,5,6,7,8,9])
    acc_prob = 1
    for i in  1:number_of_kinks
        acc_prob *= rand(add_single_kinks)(c,e)
        if acc_prob == 0
            return 0
        end
    end
    for i in  1:number_of_kinks
        acc_prob *= rand(remove_single_kinks)(c,e)
        if acc_prob == 0
            return 0
        end
    end
    return(acc_prob)
end
