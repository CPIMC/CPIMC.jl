
function move_particle(c::Configuration, e::Ensemble)

    free_orbitals = get_non_interacting_orbs_of_set(c, c.occupations)
    if length(free_orbitals) == 0
        return 1
    else
        x = rand(get_non_interacting_orbs_of_set(c, c.occupations))
    end
    oe = get_orbs_with_spin(get_non_interacting_orbs_of_set(c,setdiff!(get_sphere(x), c.occupations)), x.spin)

    #if there are no empty non interacting orbitals in neighbourhood make no change
    if length(oe) == 0
        return(1)
    end
    #@assert length(oe) > 0 "no empty orbitals in neighbourhood"

    y = rand(oe)
    @assert x != y "same Configuration proposed."

    # weight change
    dw = exp(-e.beta*(get_energy(y)-get_energy(x)))

    # change Configuration
    delete!(c.occupations,x)
    push!(c.occupations,y)

    # get orbitals for reverse update
    oe2 = setdiff!(get_sphere(y), c.occupations)

    # quotient of proposal probabilities
    dv = length(oe)/length(oe2)

    dv*dw
end

function Add_Type_B(c::Configuration, e::Ensemble)
    #samplign propability
    prop_prob = 1
    #get first tau
    tau1 = rand()*e.beta
    prop_prob *= 1/e.beta
    #get effected orbitlas (abcd)
    occs = get_occupations_at(c, tau1)
    orb_c = rand(occs)
    prop_prob *= 1/e.N
    orb_d = orb_c
    while orb_d == orb_c
        orb_d = rand(occs)
    end
    prop_prob *= 1/(e.N-1)
    opportiunisties_orb_a = get_orbs_with_spin(setdiff!(get_sphere(orb_c), occs),orb_c.spin)
    opportiunisties_orb_b = get_orbs_with_spin(setdiff!(get_sphere(orb_d), occs),orb_d.spin)
    if (length(opportiunisties_orb_a) == 0) | (length(opportiunisties_orb_b) == 0)
        return 1
    end
    orb_a = rand(opportiunisties_orb_a)
    orb_b = Orbital((orb_c.vec-orb_a.vec) + orb_d.vec,orb_d.spin)
    if (!in(orb_b,opportiunisties_orb_b) | (orb_a == orb_b))
        return 1
    end


    #there are two ways to propose the same update that we have to consider
    prop_prob = prop_prob*1/length(opportiunisties_orb_a) + prop_prob*1/length(opportiunisties_orb_b)

    #get tau2
    boarders = get_Tau_boarders(c, Set([orb_a,orb_b,orb_c,orb_d]),tau1,e)
    tau2 = rand()*(boarders[2]-boarders[1]) + boarders[1]

    if tau1 < tau2
        firsttau = tau1
        lasttau  = tau2
    else
        firsttau = tau2
        lasttau  = tau1
    end
    if lasttau > e.beta
        lasttau -= e.beta
    elseif firsttau < 0
        firsttau += e.beta
    end
    #there are 2 possibilites that will result into the same 2 taus
    prop_prob *= 2/(boarders[2]-boarders[1])

    #change configuration
    #print(T4(orb_a,orb_b,orb_c,orb_d))
    c.kinks[firsttau] = T4(orb_a,orb_b,orb_c,orb_d)
    c.kinks[lasttau] = T4(orb_c,orb_d,orb_a,orb_b)

    #Look at wether the new pair of Kinks modifies the Occupations at Tau = 0
    if firsttau > lasttau
        change_occupations(c.occupations, T4(orb_a,orb_b,orb_c,orb_d))
    end

    # quotient of proposal probabilities
    dv = length(c.kinks)/prop_prob
    # weight factor
    dw = get_abs_offdiagonal_element(e,c,T4(orb_a,orb_b,orb_c,orb_d))^2 *
            exp(-(max(tau1,tau2)-min(tau1,tau2)) * (get_energy(orb_a)
                + get_energy(orb_b) -get_energy(orb_c) -get_energy(orb_d)))
    #return quotient of poposing probabilites
    return(dv*dw)
end

function remove_Type_B(c::Configuration, e::Ensemble)
    if length(c.kinks) == 0
        return 1
    end
    Kink1 = rand(c.kinks)
    prop_prob = length(c.kinks)
    Tau_Kink2 = last(get_Tau_boarders(c, Set([last(Kink1).i, last(Kink1).j, last(Kink1).k, last(Kink1).l]),first(Kink1),e))
    if Tau_Kink2 > e.beta
        Tau_Kink2 -= e.beta
    end
    if !(Tau_Kink2 in keys(c.kinks))
        "hier breakpoint"
    end
    Kink2 = Tau_Kink2 => c.kinks[Tau_Kink2]
    #look if Kinks are type-b-connected
    ijkl = Set([last(Kink1).i, last(Kink1).j, last(Kink1).k, last(Kink1).l])
    if ijkl ==
        Set([last(Kink2).i, last(Kink2).j, last(Kink2).k, last(Kink2).l])
        #print(ijkl)
        delete!(c.kinks, first(Kink1))
        delete!(c.kinks, first(Kink2))
        #see if occupations at tau=0 are modified
        if first(Kink1) > first(Kink2)
            change_occupations(c.occupations, last(Kink2))
        end
        #calculate inverse prop_prob (see  Add_Type_B)
        boarders = get_Tau_boarders(c, ijkl ,first(Kink1),e)
        k = last(Kink1).k
        l = last(Kink1).l
        inverse_prop_prob = (1/e.N)*(1/e.N-1)*
            (1/length(get_orbs_with_spin(setdiff!(get_sphere(k), c.occupations),k.spin)) +
                1/length(get_orbs_with_spin(setdiff!(get_sphere(l), c.occupations),l.spin)))
            *(1/e.beta)*2/(boarders[2]-boarders[1])
        # quotient of proposal probabilities
        dv = inverse_prop_prob/prop_prob
        # weight factor
        dw = 1/(get_abs_offdiagonal_element(e,c,last(Kink1)))^2 *
                exp((first(Kink2)-first(Kink1)) * (get_energy(last(Kink1).i)
                    + get_energy(last(Kink1).j) -get_energy(last(Kink1).k) -get_energy(last(Kink1).l)))
        return (dv*dw)
    else
        return(1)
    end
end






function add_particle()
end

function remove_particle()
end
