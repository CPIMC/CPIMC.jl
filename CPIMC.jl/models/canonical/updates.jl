function move_particle(e::Ensemble, c::Configuration)

    x = rand(c.occupations)
    y = rand(emptyOrbs(e,c))

    delete!(c.occupations,x)
    push!(c.occupations,y)

    return exp(-e.beta*(K(y,e,c)-K(x,e,c)))
end

function reject!(c, old_c)
    copy!(c, old_c)
end

function update(e,c,updates)
    copy!(old_conf,c)
    update = rand(updates)
    if rand() < update(e,c)
        # accept
    else
        reject!(c, old_conf)
    end
end
