function Ekin(c::Configuration)
    sum(get_energy(n) for n in c.occupations)
end

function occupations(c::Configuration, emax::Int=100)
    nk = zeros(emax)

    for en in get_energy.(c.occupations)
        nk[en+1] = nk[en+1] + 1
    end
    nk
end

function particleNumber(c::Configuration)
  return length(c.occupations)
end
