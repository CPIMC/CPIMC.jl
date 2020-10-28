const ex_radius = 2 #max Radius for exitation
function move_particle(c::Configuration, e::Ensemble)
    free_orbitals = get_non_interacting_orbs_of_set(c, c.occupations)
    if length(free_orbitals) == 0
        return 1
    else
        x = rand(get_non_interacting_orbs_of_set(c, c.occupations))
    end
    oe = get_non_interacting_orbs_of_set(c,get_orbs_with_spin(setdiff!(get_sphere(x, dk = ex_radius), c.occupations), x.spin))

    #if there are no empty non interacting orbitals in neighbourhood make no change
    if length(oe) == 0
        return(1)
    end
    y = rand(oe)
    @assert x != y "same Configuration proposed."

    delta_di = get_change_diagonal_interaction(c, e, T2(y,x), img_time(0), img_time(1))
    @assert delta_di != Inf
    # weight change
    dw = exp(-(e.beta*(get_energy(y)-get_energy(x)) + delta_di))
    #println("kindiff:", e.beta*(get_energy(y)-get_energy(x)))
    #println("wdiagdiff:", delta_di,"\n")

    # change Configuration
    delete!(c.occupations,x)
    push!(c.occupations,y)

    # get orbitals for reverse update
    oe2 = get_non_interacting_orbs_of_set(c,get_orbs_with_spin(setdiff!(get_sphere(y, dk = ex_radius), c.occupations), y.spin))

    # quotient of proposal probabilities
    dv = length(oe)/length(oe2)
    @assert delta_di != Inf
    @assert (dv*dw) >= 0
    return (dv*dw)
end

function Add_Type_B(c::Configuration, e::Ensemble)
    #samplign propability
    prop_prob = 1
    #get first tau
    tau1 = img_time(rand())
    #do not allow two kinks at the same time
    while haskey(c.kinks, tau1)
        tau1 = img_time(rand())
    end
    occs = get_occupations_at(c, tau1)
    orb_c = rand(occs)
    prop_prob *= 1/e.N
    orb_d = rand(occs)
    while orb_d == orb_c
        orb_d = rand(occs)
    end
    prop_prob *= 1/(e.N-1)
    opportiunisties_orb_a = get_orbs_with_spin(setdiff!(get_sphere(orb_c, dk = ex_radius), occs),orb_c.spin)
    opportiunisties_orb_b = get_orbs_with_spin(setdiff!(get_sphere(orb_d, dk = ex_radius), occs),orb_d.spin)
    if (length(opportiunisties_orb_a) == 0) | (length(opportiunisties_orb_b) == 0)
        return 1
    end
    orb_a = rand(opportiunisties_orb_a)
    orb_b = Orbital((orb_c.vec-orb_a.vec) + orb_d.vec,orb_d.spin)
    @assert !((orb_a == orb_d) | (orb_b == orb_c))
    if (!in(orb_b,opportiunisties_orb_b) | (orb_a == orb_b))
        return 1
    end


    #we do consider kinks that differ only in the order of indizes differnt states in
    #the marcov chain although they are identical in weight and estimators
    #otherwise we would need a factor here
    prop_prob *= 1/length(opportiunisties_orb_a)

    #get tau2
    boarders = get_Tau_boarders(c, Set([orb_a,orb_b,orb_c,orb_d]),tau1)
    possible_tau2_interval = boarders[2]-boarders[1]
    if possible_tau2_interval < 0
        possible_tau2_interval += 1
    end
    tau2 = img_time(rand()*(possible_tau2_interval) + boarders[1])
    if tau2 > 1
        tau2 -= 1
    end
    #do not allow two kinks at the same time
    while haskey(c.kinks, tau2)
        tau2 = img_time(rand()*(possible_tau2_interval) + boarders[1])
        if tau2 > 1
            tau2 -= 1
        end
    end
    #See whick of the two imiginary time point is the "left boarder of the intervall"
    if tau1 > boarders[1]
        if tau2 > boarders[1]
            firsttau = min(tau1,tau2)
            lasttau = max(tau1,tau2)
        else
            firsttau = tau1
            lasttau = tau2
        end
    else
        if tau2 < boarders[1]
            firsttau = min(tau1,tau2)
            lasttau = max(tau1,tau2)
        else
            firsttau = tau2
            lasttau = tau1
        end
    end
    delta_Tau = lasttau - firsttau
    if delta_Tau < 0
        delta_Tau += 1
    end
    @assert(delta_Tau > 0)
    @assert(delta_Tau <= 1)
    #there are 2 possibilites that will result into the same 2 taus
    prop_prob *= 2.0/float(possible_tau2_interval)
    #calculate change in diagonal interaction energy
    delta_di = get_change_diagonal_interaction(c, e, T4(orb_a,orb_b,orb_c,orb_d), firsttau, lasttau)

    #change configuration
    c.kinks[firsttau] = T4(orb_a,orb_b,orb_c,orb_d)
    #shuffle index order of the second Kink
    if rand() < 0.5
        if rand() < 0.5
            #swap both
            c.kinks[lasttau] = T4(orb_d,orb_c,orb_b,orb_a)
        else
            #swap anihilators
            c.kinks[lasttau] = T4(orb_c,orb_d,orb_b,orb_a)
        end
    else
        if rand() < 0.5
            #swap creators
            c.kinks[lasttau] = T4(orb_d,orb_c,orb_a,orb_b)
        else
            #swap none
            c.kinks[lasttau] = T4(orb_c,orb_d,orb_a,orb_b)
        end
    end
    prop_prob *= 1/4

    #Look at wether the new pair of Kinks modifies the Occupations at Tau = 0
    if firsttau > lasttau
        change_occupations(c.occupations, T4(orb_a,orb_b,orb_c,orb_d))
    end

    # quotient of proposal probabilities
    dv = (1/length(c.kinks))/prop_prob

    # weight factor
    dw = ((e.beta)^2) *
            get_abs_offdiagonal_element(e,c,T4(orb_a,orb_b,orb_c,orb_d))^2 *
            exp(-((delta_Tau)*e.beta * (get_energy(orb_a) +
                 get_energy(orb_b) - get_energy(orb_c) - get_energy(orb_d)) + delta_di))
    @assert (dv*dw) >= 0
    return(dv*dw)
end

function remove_Type_B(c::Configuration, e::Ensemble)

    if length(c.kinks) == 0
        return 1
    end
    Kink1 = rand(c.kinks)
    # if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it
    if dot(last(Kink1).i.vec-last(Kink1).k.vec, last(Kink1).i.vec-last(Kink1).k.vec) > ex_radius^2
        return 1
    end
    prop_prob = 1/length(c.kinks)
    Tau_Kink2 = last(get_Tau_boarders(c, Set([last(Kink1).i, last(Kink1).j, last(Kink1).k, last(Kink1).l]),first(Kink1)))
    Kink2 = Tau_Kink2 => c.kinks[Tau_Kink2]
    #look if Kinks are type-b-connected
    ijkl = Set([last(Kink1).i, last(Kink1).j, last(Kink1).k, last(Kink1).l])
    if ijkl ==
        Set([last(Kink2).i, last(Kink2).j, last(Kink2).k, last(Kink2).l])
        delete!(c.kinks, first(Kink1))
        delete!(c.kinks, first(Kink2))
        #see if occupations at tau=0 are modified
        if first(Kink1) > first(Kink2)
            change_occupations(c.occupations, last(Kink2))
            delta_Tau = first(Kink2)-first(Kink1) + 1
        else
            delta_Tau = first(Kink2)-first(Kink1)
        end
        @assert(delta_Tau > 0)
        @assert(delta_Tau <= 1)
        #calculate inverse prop_prob (see  Add_Type_B)
        boarders = get_Tau_boarders(c, ijkl ,first(Kink1))
        possible_tau2_interval = boarders[2]-boarders[1]
        if possible_tau2_interval < 0
            possible_tau2_interval = 1 + possible_tau2_interval
        end
        occs = get_occupations_at(c, first(Kink1))
        orb_a = last(Kink1).i
        orb_b = last(Kink1).j
        orb_c = last(Kink1).k
        orb_d = last(Kink1).l
        inverse_prop_prob = (1/e.N)*(1/(e.N-1)) *
            (1/length(get_orbs_with_spin(setdiff!(get_sphere(orb_c, dk = ex_radius), occs),orb_c.spin))) *
             2.0/float(possible_tau2_interval) * (1/4)
        # quotient of proposal probabilities
        dv = inverse_prop_prob/prop_prob

        #calculate change in diagonal interaction energy
        delta_di = get_change_diagonal_interaction(c, e, last(Kink1), first(Kink1), first(Kink2))


        # weight factor
        dw = (1.0/(e.beta)^2) *
            (1.0/(get_abs_offdiagonal_element(e,c,last(Kink1)))^2) *
                exp((delta_Tau)*e.beta * (get_energy(orb_a) +
                     get_energy(orb_b) - get_energy(orb_c) - get_energy(orb_d)) + delta_di)
        @assert (dv*dw) >= 0
        return (dv*dw)
    else
        return(1)
    end
end


function change_type_B(c::Configuration, e::Ensemble)
    if length(c.kinks) == 0
        return 1
    end
    #print(length(c.kinks),"\n")
    Kink1 = rand(c.kinks)
    Tau_Kink2 = last(get_Tau_boarders(c, Set([last(Kink1).i, last(Kink1).j, last(Kink1).k, last(Kink1).l]),first(Kink1)))
    Kink2 = Tau_Kink2 => c.kinks[Tau_Kink2]
    #look if Kinks are type-b-connected
    ijkl = Set([last(Kink1).i, last(Kink1).j, last(Kink1).k, last(Kink1).l])
    if ijkl != Set([last(Kink2).i, last(Kink2).j, last(Kink2).k, last(Kink2).l])
        return(1)
    end
    occs = get_occupations_at(c, first(Kink1))

    opportunities = get_non_interacting_orbs_of_set_in_interval(
                        c,get_orbs_with_spin(
                            setdiff!(
                                get_sphere(last(Kink1).i, dk = ex_radius
                                ), occs
                            ), last(Kink1).i.spin
                        ),first(Kink1),first(Kink2)
                    )
    delete!(opportunities, last(Kink1).k)
    delete!(opportunities, last(Kink1).l)
    if length(opportunities) == 0
        return 1
    end
    new_orb_i = rand(opportunities)
    new_orb_j = Orbital(last(Kink1).j.vec + last(Kink1).i.vec - new_orb_i.vec, last(Kink1).j.spin)
    if new_orb_i == new_orb_j
        return 1
    end
    if (!is_non_interacting_in_interval(c,new_orb_j,first(Kink1),first(Kink2))) |
        in(new_orb_j,occs)
        return 1
    else
        #calculate change in diagonal interaction energy
        delta_di = get_change_diagonal_interaction(c, e, T4(new_orb_i, new_orb_j, last(Kink1).i, last(Kink1).j), first(Kink1), first(Kink2))

        #change occupations
        if first(Kink1) > first(Kink2)
            change_occupations(c.occupations, T4(new_orb_i, new_orb_j, last(Kink1).i, last(Kink1).j))
            delta_Tau = first(Kink2)-first(Kink1) + 1
        else
            delta_Tau = first(Kink2)-first(Kink1)
        end


        dw = exp(-(e.beta*delta_Tau*(get_energy(new_orb_i) + get_energy(new_orb_j) -
                                    get_energy(last(Kink1).i) - get_energy(last(Kink1).j)) + delta_di)) *
            (get_abs_offdiagonal_element(e,c,T4(new_orb_i, new_orb_j, last(Kink1).k, last(Kink1).l))/
                        get_abs_offdiagonal_element(e,c,last(Kink1)))^2

        #change Kinks
        c.kinks[first(Kink1)] = T4(new_orb_i, new_orb_j, last(Kink1).k, last(Kink1).l)
        c.kinks[first(Kink2)] = T4(last(Kink2).i, last(Kink2).j, new_orb_i, new_orb_j,)
        change_occupations(occs, T4(new_orb_i, new_orb_j, last(Kink1).i, last(Kink1).j))
        opportunites_reverse = get_non_interacting_orbs_of_set_in_interval(
                                    c,get_orbs_with_spin(
                                        setdiff!(
                                            get_sphere(new_orb_i, dk = ex_radius
                                            ), occs
                                        ), last(Kink1).i.spin
                                    ),first(Kink1),first(Kink2)
                                )
        delete!(opportunites_reverse, last(Kink1).k)
        delete!(opportunites_reverse, last(Kink1).l)
        @assert (dw * length(opportunities)/length(opportunites_reverse)) >= 0
        return (dw * length(opportunities)/length(opportunites_reverse))
    end
end

function shuffle_indixes(c::Configuration, e::Ensemble)
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



function add_particle()
end

function remove_particle()
end
