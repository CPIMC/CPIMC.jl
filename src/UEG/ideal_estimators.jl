#TODO Delete this File?
function Ekin(e::Ensemble, c::Configuration) :: UInt
    sum(energy(n) for n in c.occupations)
end

function occupations(e::Ensemble, c::Configuration, emax::Int=100) :: Array{UInt,1}
    nk = zeros(UInt, emax)

    ens = energy.(c.occupations)

    for ε in ens[ens .< emax]
        nk[ε+1] = nk[ε+1] + 1
    end
    nk
end

function particleNumber(e::Ensemble, c::Configuration)
  return length(c.occupations)
end

