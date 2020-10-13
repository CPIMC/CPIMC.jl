function Ekin(c::Configuration)
    sum(get_energy(n) for n in c.occupations)
end

function occVec(c::Configuration)
    get_energy.(c.occupations)
end

function occupations(c::Configuration, emax::Int=99)
    nk = zeros(emax)
    for en in 0:emax
        if in(en, get_energy.(c.occupations))
            nk[en+1] = 1
        end
    end
    nk
end

function particleNumber(c::Configuration)
  return length(c.occupations)
end
