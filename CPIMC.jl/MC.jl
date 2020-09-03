function sweep(steps::Int, sampleEvery::Int, updates, measurements, e::Ensemble, c::Configuration)

    i = 0

    while i < steps

        update(e,c,updates)

        if i % sampleEvery == 0
            for (stat,obs) in measurements
                if typeof(stat) == Group
                    fit!(stat, eachrow(obs(e,c)))
                else
                    fit!(stat, obs(e,c))
                end
            end
        end

    i += 1
  end

  return zip(updatesAccepted,updatesProposed)
end
