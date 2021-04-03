
"""Returns True if left_kink and right_kink are entangled in a Type-D way
This does not check wether the two kinks are neighbouring"""
function is_type_D(left_kink::T4, right_kink::T4)
  if !(Set([left_kink.i, left_kink.j]) == Set([right_kink.k, right_kink.l])) &
        (Set([left_kink.k, left_kink.l]) == Set([right_kink.i, right_kink.j]))
    return(true)
  else
    return(false)
  end
end

"""Return a Tuple of 2 imaginaty times of 'neighbouring' Kinks that are Type-D-Entangeld AND removable.
'neighbouring' refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-D-entanglement is oriented
#to the left of the first τ."""
function get_left_type_D_removable_pairs(c::Configuration)
  pairs_left = Set{Tuple{Fixed{Int64,60},Fixed{Int64,60}}}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = τ_borders(c, kink_orb_set ,τ)
    if is_type_D(c.kinks[τ_left], kink)
      if dot(kink.i.vec-kink.k.vec,kink.i.vec-kink.k.vec) <= (ex_radius^2)
        if kink.i.spin == kink.k.spin
          push!(pairs_left, (τ, τ_left))
        end
      end
    end
  end
  return pairs_left
end

"""Return a Tuple of 2 imaginary times of 'neighbouring' Kinks that are Type-D-Entangeld AND removable.
'neighbouring' refers to that only Tuples of Kinks that are the closest Kink to act on an orbital of the
other kink in the corresponding direktion are looked at.
The Tuples are always arranged in a way that the Kink who gets neighboured by
the opther stands first.(vice versa does not have to be the case)
The Set consists of the pairs where the Type-D-entanglement is oriented
#to the right of the first τ."""
function get_right_type_D_removable_pairs(c::Configuration)
  pairs_right = Set{Tuple{Fixed{Int64,60},Fixed{Int64,60}}}()
  for (τ,kink) in c.kinks
    kink_orb_set = Set([kink.i, kink.j, kink.k, kink.l])
    τ_left,τ_right = τ_borders(c, kink_orb_set ,τ)
    if is_type_D(kink, c.kinks[τ_right])
      if dot(kink.i.vec-kink.k.vec,kink.i.vec-kink.k.vec) <= (ex_radius^2)
        if kink.i.spin == kink.k.spin
          push!(pairs_right, (τ, τ_right))
        end
      end
    end
  end
  return pairs_right
end

function possible_new_orb1_D(occs, exite_orb1, exite_orb2, old_orb1, old_orb2)
    opprtunities = filter(new_orb_1 -> in(PlaneWave(old_orb1.vec + old_orb2.vec - new_orb_1.vec, exite_orb2.spin),occs) && (new_orb_1 != PlaneWave(old_orb1.vec + old_orb2.vec - new_orb_1.vec, exite_orb2.spin)),
                    intersect!(sphere_with_same_spin(exite_orb1, dk = ex_radius), occs))

    return setdiff!(opprtunities, Set([exite_orb1, exite_orb2, old_orb1, old_orb2]))
end

function add_type_D(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64,Step}
    prop_prob = 1.0
    if isempty(c.kinks)
        return 1.0, Step()
    end
    old_kink = rand(c.kinks)
    prop_prob *= 1.0/length(c.kinks)
    occs = occupations(c, first(old_kink))
    prop_prob *= 0.5 #left or right
    if rand() >= 0.5
        #add kink left
        opportunities_new_orb1 = possible_new_orb1_D(occs, last(old_kink).i, last(old_kink).j, last(old_kink).k, last(old_kink).l)
        if isempty(opportunities_new_orb1)
            return 1.0, Step()
        end
        new_orb1 = rand(opportunities_new_orb1)
        prop_prob *= 1.0/length(opportunities_new_orb1)
        new_orb2 = PlaneWave(last(old_kink).k.vec + (last(old_kink).l.vec - new_orb1.vec), last(old_kink).j.spin)
        @assert(in(new_orb2, occs) & (new_orb1 != new_orb2))
        τ_Intervall = first(old_kink) - first(τ_borders(c, Set([last(old_kink).i, last(old_kink).j, new_orb1, new_orb2]), first(old_kink)))

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

        delta_di = ΔWdiag_element(m, e, c, last(old_kink).i,last(old_kink).j, new_orb1, new_orb2, τ_new_kink, first(old_kink))

        #change_Configuration
        #see if c.occupations change
        if τ_new_kink > first(old_kink)  #new kink was added left of old kink
            drop_orbs = Set([new_orb1, new_orb2])
            add_orbs = Set([last(old_kink).i,last(old_kink).j])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        prop_prob *= 0.5#shuffle creators of the changed kink TODO: is this necessary for ergodicy ?
        add_kinks =(
                    τ_new_kink => T4(last(old_kink).i, last(old_kink).j,new_orb1, new_orb2),
                    first(old_kink) => shuffle_creators(T4(new_orb2, new_orb1, last(old_kink).k, last(old_kink).l))
                    )
        @assert(is_type_D(last(first(add_kinks)), last(last(add_kinks))))
        drop_kinks = (old_kink,)

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))

        #calculate weight differance
        dw_off_diag = abs(Woffdiag_element(m, e, apply_step(c,Δ).kinks[τ_new_kink])) * abs(Woffdiag_element(m, e, apply_step(c,Δ).kinks[first(old_kink)])) /
                                                abs(Woffdiag_element(m, e, last(old_kink)))
        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(energy(m, last(old_kink).i) + energy(m, last(old_kink).j)-
                                                         energy(m, new_orb1) - energy(m, new_orb2)) + e.β * delta_di))

        inverse_prop_prob = (1.0/length(get_right_type_D_removable_pairs(apply_step(c,Δ)))) * 0.5
    else
        #add kink right
        opportunities_new_orb1 = possible_new_orb1_D(occs, last(old_kink).k, last(old_kink).l, last(old_kink).i, last(old_kink).j)
        if isempty(opportunities_new_orb1)
            return 1.0, Step()
        end
        new_orb1 = rand(opportunities_new_orb1)
        prop_prob *= 1.0/length(opportunities_new_orb1)
        new_orb2 = PlaneWave(last(old_kink).i.vec + (last(old_kink).j.vec - new_orb1.vec), last(old_kink).l.spin)

        @assert(in(new_orb2, occs) & (new_orb1 != new_orb2))
        τ_Intervall = last(τ_borders(c, Set([
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

        prop_prob *= 1.0/Float64(τ_Intervall)

        delta_di = ΔWdiag_element(m, e, c, last(old_kink).k, last(old_kink).l, new_orb1, new_orb2, first(old_kink), τ_new_kink)

        #change_Configuration
        #see if c.occupations change
        if τ_new_kink < first(old_kink)  #new kink was added right of old kink
            drop_orbs = Set([new_orb1,new_orb2])
            add_orbs = Set([last(old_kink).k, last(old_kink).l])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        prop_prob *= 0.5#shuffle annihilators (!) of the changed kink TODO: is this necessary for ergodicy ?
        add_kinks = (
                    τ_new_kink => T4(new_orb1, new_orb2, last(old_kink).k, last(old_kink).l),
                    first(old_kink) => shuffle_annihilators(T4(last(old_kink).i, last(old_kink).j, new_orb2, new_orb1))
                    )
        drop_kinks = (old_kink,)
        @assert(is_type_D(last(last(add_kinks)), last(first(add_kinks))))
        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))

        dw_off_diag = abs(Woffdiag_element(m, e, apply_step(c,Δ).kinks[τ_new_kink])) * abs(Woffdiag_element(m, e, apply_step(c,Δ).kinks[first(old_kink)])) /
                                                abs(Woffdiag_element(m, e, last(old_kink)))

        dw = e.β * dw_off_diag* exp(-(e.β * delta_τ*(energy(m,last(old_kink).k) + energy(m,last(old_kink).l) -
                                                         energy(m,new_orb1) - energy(m,new_orb2)) + e.β * delta_di))

        inverse_prop_prob = (1.0/length(get_left_type_D_removable_pairs(apply_step(c,Δ)))) * 0.5
    end

    @assert(delta_τ > 0 )
    @assert(!isinf((inverse_prop_prob/prop_prob)*dw))
    return ((inverse_prop_prob/prop_prob)*dw), Δ
end

function remove_type_D(m::Model, e::Ensemble, c::Configuration) :: Tuple{Float64,Step}
    prop_prob = 0.5
    if rand() > 0.5
        #removed kink left of changed kink
        opportunities = get_right_type_D_removable_pairs(c)
        if isempty(opportunities)
            return 1.0, Step()
        end
        removed_kink_τ, changed_kink_τ = rand(opportunities)
        prop_prob *= 1.0/length(opportunities)

        @assert (c.kinks[removed_kink_τ].i.spin == c.kinks[removed_kink_τ].k.spin)
        @assert (dot(c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec,c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec) <= (ex_radius^2)) "if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it"



        #safe those for later
        removed_orb1 = c.kinks[removed_kink_τ].k
        removed_orb2 = c.kinks[removed_kink_τ].l

        #change configuration
        @assert removed_orb1 != c.kinks[changed_kink_τ].k
        add_kinks = (changed_kink_τ => T4(c.kinks[removed_kink_τ].i,c.kinks[removed_kink_τ].j,c.kinks[changed_kink_τ].k,c.kinks[changed_kink_τ].l),)

        #see if c.occupations change
        if removed_kink_τ > changed_kink_τ
            # change_occupations(c.occupations, T4(removed_orb1,removed_orb2, c.kinks[changed_kink_τ].i,c.kinks[changed_kink_τ].j))
            drop_orbs = Set([c.kinks[removed_kink_τ].i,c.kinks[removed_kink_τ].j])
            add_orbs = Set([removed_orb1,removed_orb2])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end


        # delete!(c.kinks, removed_kink_τ)
        drop_kinks = (removed_kink_τ => c.kinks[removed_kink_τ], changed_kink_τ => c.kinks[changed_kink_τ], )

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))

        #calculate reverse_prop_prob
        occs = occupations(apply_step(c,Δ), changed_kink_τ)
        opportunities_reverse_new_orb1 = possible_new_orb1_D(occs, apply_step(c,Δ).kinks[changed_kink_τ].i, apply_step(c,Δ).kinks[changed_kink_τ].j, apply_step(c,Δ).kinks[changed_kink_τ].k, apply_step(c,Δ).kinks[changed_kink_τ].l)
        @assert(in(removed_orb1,opportunities_reverse_new_orb1))
        τ_Intervall = changed_kink_τ - first(τ_borders(apply_step(c,Δ), Set([removed_orb1, removed_orb2,
                         apply_step(c,Δ).kinks[changed_kink_τ].i, apply_step(c,Δ).kinks[changed_kink_τ].j]),changed_kink_τ))


        if τ_Intervall < 0
            τ_Intervall +=1
        end

        inverse_prop_prob = (0.5/length(apply_step(c,Δ).kinks)) * (1.0/length(opportunities_reverse_new_orb1)) *
                                 (1.0/Float64(τ_Intervall)) * (1/2)# TODO: (1/2) = 0.5

        #calculate weight change
        delta_di = ΔWdiag_element(m, e, apply_step(c,Δ), apply_step(c,Δ).kinks[changed_kink_τ].i,apply_step(c,Δ).kinks[changed_kink_τ].j, removed_orb1,removed_orb2, removed_kink_τ, changed_kink_τ)

        dw_off_diag = abs(Woffdiag_element(m, e, T4(removed_orb1, removed_orb2, apply_step(c,Δ).kinks[changed_kink_τ].k, apply_step(c,Δ).kinks[changed_kink_τ].l))) *
                        abs(Woffdiag_element(m, e,T4(apply_step(c,Δ).kinks[changed_kink_τ].i, apply_step(c,Δ).kinks[changed_kink_τ].j, removed_orb1, removed_orb2))) /
                            abs(Woffdiag_element(m, e,apply_step(c,Δ).kinks[changed_kink_τ]))

        delta_τ = Float64(changed_kink_τ - removed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1.0/e.β)* (1.0/dw_off_diag) * exp(e.β * delta_τ * (energy(m,apply_step(c,Δ).kinks[changed_kink_τ].i) + energy(m,apply_step(c,Δ).kinks[changed_kink_τ].j) -
                                                                    energy(m,removed_orb1) - energy(m,removed_orb2)) + e.β * delta_di)

    else
        #removed kink right of changed kink
        opportunities = get_left_type_D_removable_pairs(c)
        if isempty(opportunities)
            return 1.0, Step()
        end
        removed_kink_τ, changed_kink_τ = rand(opportunities)
        prop_prob *= 1.0/length(opportunities)


        @assert (c.kinks[removed_kink_τ].i.spin == c.kinks[removed_kink_τ].k.spin)
        @assert (dot(c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec,
                    c.kinks[removed_kink_τ].i.vec-c.kinks[removed_kink_τ].k.vec) <= (ex_radius^2)) "if the difference between i and k is larger then ex_radius we can not create the kink and therefore also can't delete it"


        #safe thoose for later
        removed_orb1 = c.kinks[removed_kink_τ].i
        removed_orb2 = c.kinks[removed_kink_τ].j


        #change configuration
        add_kinks = (changed_kink_τ => T4(c.kinks[changed_kink_τ].i,c.kinks[changed_kink_τ].j,c.kinks[removed_kink_τ].k,c.kinks[removed_kink_τ].l),)

        #see if c.occupations change
        if removed_kink_τ < changed_kink_τ  #new kink was added right of old kink
            drop_orbs = Set([c.kinks[removed_kink_τ].k,c.kinks[removed_kink_τ].l])
            add_orbs = Set([removed_orb1, removed_orb2])
        else
            drop_orbs = Set{basis(c)}()
            add_orbs = Set{basis(c)}()
        end

        # delete!(c.kinks, removed_kink_τ)
        drop_kinks = (removed_kink_τ => c.kinks[removed_kink_τ], changed_kink_τ => c.kinks[changed_kink_τ])

        # MC Step generated by this update
        Δ = Step(Configuration(drop_orbs, drop_kinks...), Configuration(add_orbs, add_kinks...))

        #calculate reverse_prop_prob
        occs = occupations(apply_step(c,Δ), changed_kink_τ)
        opportunities_reverse_new_orb1 = possible_new_orb1_D(occs, apply_step(c,Δ).kinks[changed_kink_τ].k, apply_step(c,Δ).kinks[changed_kink_τ].l, apply_step(c,Δ).kinks[changed_kink_τ].i, apply_step(c,Δ).kinks[changed_kink_τ].j)
        @assert(in(removed_orb1,opportunities_reverse_new_orb1))
        τ_Intervall =  last(τ_borders(apply_step(c,Δ), Set([apply_step(c,Δ).kinks[changed_kink_τ].k, apply_step(c,Δ).kinks[changed_kink_τ].l,
                                                    removed_orb1, removed_orb2]),changed_kink_τ)) - changed_kink_τ

        if τ_Intervall < 0
            τ_Intervall +=1
        end
        inverse_prop_prob = (0.5/length(apply_step(c,Δ).kinks)) * (1.0/length(opportunities_reverse_new_orb1)) *
                                 (1.0/Float64(τ_Intervall)) * (1/2)# TODO: (1/2) = 0.5

        #calculate weight change
        delta_di = ΔWdiag_element(m, e, apply_step(c,Δ), apply_step(c,Δ).kinks[changed_kink_τ].k,apply_step(c,Δ).kinks[changed_kink_τ].l, removed_orb1,removed_orb2, changed_kink_τ, removed_kink_τ)

        dw_off_diag = abs(Woffdiag_element(m, e, T4(removed_orb1, removed_orb2,  apply_step(c,Δ).kinks[changed_kink_τ].k, apply_step(c,Δ).kinks[changed_kink_τ].l))) *
                        abs(Woffdiag_element(m, e, T4(apply_step(c,Δ).kinks[changed_kink_τ].i, apply_step(c,Δ).kinks[changed_kink_τ].j, removed_orb1, removed_orb2))) /
                            abs(Woffdiag_element(m, e, apply_step(c,Δ).kinks[changed_kink_τ]))

        delta_τ = Float64(removed_kink_τ - changed_kink_τ)
        if delta_τ < 0
            delta_τ +=1
        end
        dw = (1.0/e.β)*(1.0/dw_off_diag) * exp(e.β * delta_τ*(energy(m,apply_step(c,Δ).kinks[changed_kink_τ].k) + energy(m, apply_step(c,Δ).kinks[changed_kink_τ].l) -
                                                            energy(m, removed_orb1) - energy(m, removed_orb2)) + e.β * delta_di)

    end

    @assert(dw != Inf)
    @assert(delta_τ > 0 )
    @assert(!isinf((inverse_prop_prob/prop_prob) * dw))

    return ((inverse_prop_prob/prop_prob) * dw), Δ
end
