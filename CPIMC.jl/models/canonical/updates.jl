function move_particle(e::Ensemble, c::Configuration)

    x = rand(c.occupations)
    y = rand(emptyOrbs(e,c))

    delete!(c.occupations,x)
    push!(c.occupations,y)

    return exp(-e.beta*(get_energy(y,orblist)-get_energy(x,orblist)))
end

function update(e,c,updates)
    @assert !iszero(length(updates))
    old_conf = Configuration(copy(c.occupations))#Funktion in klasse einbauen?
    up = rand(updates)
    dp = up(e,c)
    if rand() < dp
        # accept
    else
        # reject
        ## TODO: is this a new object ???
        c.occupations = old_conf.occupations
    end
end
