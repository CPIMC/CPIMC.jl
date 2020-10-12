function Ekin(e::Ensemble, c::Configuration)
    sum(get_energy(n) for n in c.occupations)
end

function occVec(e::Ensemble, c::Configuration)
    get_energy.(c.occupations)
end

function particleNumber(e::Ensemble, c::Configuration)
  return length(c.occupations)
end
