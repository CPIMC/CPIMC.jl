function sweep(steps::Int, sampleEvery::Int, updates, measurements, e::Ensemble, c::Configuration)
  nUpdates = length(updates)

  @assert nUpdates > 0

  updatesProposed = zeros(Int, nUpdates)
  updatesAccepted = zeros(Int, nUpdates)

  i = 0

  while i < steps
    k = rand(1:nUpdates)

    result = updates[k](e,c)

    updatesProposed[k] += 1
    
    if result == :drop
      continue
    end

    if result == :accept
      updatesAccepted[k] += 1
    end

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

