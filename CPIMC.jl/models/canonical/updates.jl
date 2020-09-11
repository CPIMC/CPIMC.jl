function move_particle(e::Ensemble, c::Configuration, orblist)

    x = rand(c.occupations)
    y = rand(emptyOrbs(e,c))

    delete!(c.occupations,x)
    push!(c.occupations,y)

    return exp(-e.beta*(get_energy(y,orblist)-get_energy(x,orblist)))
end

function update(e,c,updates,orblist)
    @assert !iszero(length(updates))
    old_conf = c
    up = rand(updates)
    dp = up(e,c,orblist)
    if rand() < dp
        # accept
    else
        # reject
        ## TODO: is this a new object ???
        c = old_conf
    end
end
