function move_particle(c::Configuration)
    x,y = proposeOrbs(c)

    # weight change
    dw = exp(-e.beta*(get_energy(y,orblist)-get_energy(x,orblist)))

    if rand() < dw
        # accept
        delete!(c.occupations,x)
        push!(c.occupations,y))
    else
        # reject
        ## TODO: is this a new object ???
        # c.occupations = old_conf.occupations
    end
end

function update(c,updates)
    @assert !iszero(length(updates))
    up = rand(updates)
    up(c)
end
