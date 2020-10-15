function Ekin(c::Configuration)
    sum(get_energy(n) for n in c.occupations)
end

function occupations(c::Configuration, emax::Int=100)
    nk = zeros(emax)

    ens = get_energy.(c.occupations)
    ens = ens[ens .< emax]

    for ε in ens
        nk[ε+1] = nk[ε+1] + 1
    end
    nk
end

function particleNumber(c::Configuration)
  return length(c.occupations)
end
