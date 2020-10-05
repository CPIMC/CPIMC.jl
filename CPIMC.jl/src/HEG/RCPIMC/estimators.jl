function Ekin(e::Ensemble, c::Configuration)
    if length(c.kinks) == 0
        return(sum(get_energy(n) for n in c.occupations))
    end
    occs = copy(c.occupations)
    E_kin = 0
    #betrachte das komplette intervall zwischen letztem und erstem Kink als
    #erstes intervall
    old_Tau = first(last(c.kinks)) - e.beta
    for (Tau,kink) in c.kinks
        E_kin += sum(get_energy(n) for n in c.occupations) * (Tau - old_Tau)
        old_Tau = Tau
        change_occupations(occs, kink)
    end
    return(E_kin * 1/e.beta)
end

function W_off_diag(e::Ensemble, c::Configuration)
    return(-(length(c.kinks)/e.beta))
end

function W_Diag(e::Ensemble, c::Configuration)
    W_diag = 0
    if length(c.kinks) == 0

        for occ1 in c.occupation
            redundant = true
            for occ2 in c.occupation
                if !redundant
                    W_diag += (lamda(e.N,e.rs)/2) * 1/dot((occ1.vec-occ2.vec),(occ1.vec-occ2.vec))
                end
                if occ1 == occ2
                    redundant == false
                end
            end
        end
    end
end

function E_pot(e::Ensemble, c::Configuration)
    return(W_Diag(e, c) + W_off_diag(e,c))
end

############Noch nicht fürs WW-System
function occVec(e::Ensemble, c::Configuration)
    get_energy.(c.occupations)
end

function particleNumber(c::Configuration)
  return length(c.occupations)
end
