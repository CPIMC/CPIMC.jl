function move_particle(e::Ensemble, c::Configuration, orblist)

    x = rand(c.occupations)
    y = rand(emptyOrbs(e,c))

    delete!(c.occupations,x)
    push!(c.occupations,y)

    return exp(-e.beta*(Ekin(y,e,c,orblist)-Ekin(x,e,c,orblist)))
end

function update(e,c,updates)
    @assert !iszero(length(updates))
    old_conf = c
    up = rand(updates)
    dp = up(e,c)
    if rand() < dp
        # accept
    else
        # reject
        ## TODO: is this a new object ???
        c = old_conf
    end
end
