function move_particle(e::Ensemble, c::Configuration, orblist)

    x = rand(c.occupations)
    y = rand(emptyOrbs(e,c))

    delete!(c.occupations,x)
    push!(c.occupations,y)

    return 1 #exp(-e.beta*(Ekin(y,e,c,orblist)-Ekin(x,e,c,orblist)))
end

function update(e,c,updates,orblist, max_update_length)
    @assert !iszero(length(updates))
    update_length =  rand(1:max_update_length)
    old_conf = c
    sampling_prob_quod = 1
    for step_index in update_length
        up = rand(updates) 
        sampling_prob_quod *= up(e,c,orblist)
    end
    acc_prob = sampling_prob_quod * exp(-e.beta*(Ekin(e,old_conf,orblist)-Ekin(e,c,orblist)))

    if rand() < acc_prob
        # accept
    else
        # reject
        ## TODO: is this a new object ???
        c = old_conf

    end
end
