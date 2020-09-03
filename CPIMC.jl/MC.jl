function sweep(steps::Int, sampleEvery::Int, updates, measurements, e::Ensemble, c::Configuration, orblist)

    i = 0

    while i < steps

        update(e,c,updates,orblist)

        if i % sampleEvery == 0
            for (stat,obs) in measurements
                if typeof(stat) == Group
                    fit!(stat, eachrow(obs(e,c,orblist)))
                else
                    fit!(stat, obs(e,c,orblist))
                end
            end
        end

    i += 1
  end

end
