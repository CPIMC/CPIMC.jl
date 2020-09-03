function add_particle(e::Ensemble, c::Configuration)
  @assert e.cutoff > c.N + 1

  x = rand(emptyOrbs(e,c))

  if rand() < exp(-e.beta*(K(x,e,c)-e.mu))
    push!(c.occupations,x)
    c.N += 1
    return :accept
  else
    return :reject
  end
end

function remove_particle(e::Ensemble, c::Configuration)
  if c.N == 0
    return :drop
  end
 
  x = rand(c.occupations)

  if rand() < exp(-e.beta*(-K(x,e,c)+e.mu))
    delete!(c.occupations,x)
    c.N -= 1
    return :accept
  else
    return :reject
  end
end

function move_particle(e::Ensemble, c::Configuration)
  if c.N == 0
    return :drop
  end

  x = rand(c.occupations)
  y = rand(emptyOrbs(e,c))

  if rand() < exp(-e.beta*(K(y,e,c)-K(x,e,c)))
    delete!(c.occupations,x)
    push!(c.occupations,y)

    return :accept
  else
    return :reject
  end
end
