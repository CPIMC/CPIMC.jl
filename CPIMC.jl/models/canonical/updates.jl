function move_particle(e::Ensemble, c::Configuration)

    x = rand(c.occupations)
    y = rand(emptyOrbs(e,c))

    delete!(c.occupations,x)
    push!(c.occupations,y)

    return exp(-e.beta*(K(y,e,c)-K(x,e,c)))
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
