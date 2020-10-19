function move_particle!(c::Configuration, e::Ensemble) :: Tuple{Float64, Step}
    x = rand(c.occupations)
    oe = setdiff!(get_sphere(x), c.occupations)
    if isempty(oe)
        return 1, Step()
    else
        y = rand(oe)
        @assert x != y "same Configuration proposed."

        # get orbitals for reverse update
        oe2 = setdiff!(get_sphere(y), c.occupations)

        # quotient of proposal probabilities
        δv = length(oe)/length(oe2)

        return δv * weight(y, x, e), Step(y,x)
    end
end

function add_particle!()
end

function remove_particle!()
end
