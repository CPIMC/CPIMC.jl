function move_particle(e::Ensemble, c::Configuration)

    x = rand(c.occupations)
    y = rand(emptyOrbs(e,c))

    delete!(c.occupations,x)
    push!(c.occupations,y)

    return exp(-e.beta*(K(y,e,c)-K(x,e,c)))
end

# TODO : Avoid fields with abstract type
struct Updates
    Set([move_particle])
end

function reject!(c, old_c)
    copy!(c, old_c)
end

function update(e,c)
    copy!(old_conf,c)
    update = rand(Updates)
    if rand() < update(e,c)
        # accept
    else
        reject!(c, old_conf)
    end
end
