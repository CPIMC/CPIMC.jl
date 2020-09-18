function sweep(steps::Int, sampleEvery::Int, throwAway::Int, updates, measurements, e::Ensemble, c::Configuration)

    for i in 1:throwAway
        update!(c,e,updates)
    end

    i = 0

    while i < steps
        update!(c,e,updates)

        if i % sampleEvery == 0
            for (key,(stat,obs)) in measurements
                if typeof(stat) == Group
                    fit!(stat, eachrow(obs(c)))
                else
                    fit!(stat, obs(c))
                end
            end
        end

        i += 1
    end
end
