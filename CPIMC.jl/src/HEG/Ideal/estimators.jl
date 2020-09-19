function Ekin(c::Configuration)
    sum(get_energy(n) for n in c.occupations)
end

function occVec(c::Configuration)
    get_energy.(c.occupations)
end

function particleNumber(c::Configuration)
  return length(c.occupations)
end
