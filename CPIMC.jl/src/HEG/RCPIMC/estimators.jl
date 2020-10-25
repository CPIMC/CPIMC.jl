function Ekin(e::Ensemble, c::Configuration)
    if length(c.kinks) == 0
        return(sum(get_energy(n) for n in c.occupations))
    end
    occs = copy(c.occupations)
    E_kin = 0
    #betrachte das komplette intervall zwischen letztem und erstem Kink als
    #erstes intervall
    old_Tau = first(last(c.kinks)) - 1
    for (Tau,kink) in c.kinks
        E_kin += sum(get_energy(n) for n in occs) * float(Tau - old_Tau)
        old_Tau = Tau
        change_occupations(occs, kink)
    end
    return(E_kin)
end

function W_off_diag(e::Ensemble, c::Configuration)
    return(-(length(c.kinks)/e.beta))
end

function W_diag(e::Ensemble, c::Configuration)
    W_diag = 0
    if length(c.kinks) == 0
        for occ1 in c.occupations
            redundant = true
            for occ2 in c.occupations
                if !redundant
                    # W_diag += (lambda(e.N,e.rs)/2) * 1/dot((occ1.vec-occ2.vec),(occ1.vec-occ2.vec))
                    W_diag += 0.5 * lambda(e.N,e.rs) / dot((occ1.vec-occ2.vec),(occ1.vec-occ2.vec))
                end
                if occ1 == occ2
                    redundant = false
                end
            end
        end
    else
        occs = copy(c.occupations)
        old_Tau = first(last(c.kinks)) - 1
        for (Tau,kink) in c.kinks
            for occ1 in occs
                redundant = true
                for occ2 in occs
                    if !redundant
                        # W_diag += (lambda(e.N,e.rs)/2) * 1/dot((occ1.vec-occ2.vec),(occ1.vec-occ2.vec)) * (Tau-old_Tau)
                        W_diag += 0.5 * lambda(e.N,e.rs) / dot((occ1.vec-occ2.vec),(occ1.vec-occ2.vec)) * (Tau-old_Tau)
                    end
                    if occ1 == occ2
                        redundant = false
                    end
                end
            end
            old_Tau = Tau
            change_occupations(occs, kink)
        end
    end
    W_diag *= -1
    return(W_diag)
end

function Epot(e::Ensemble, c::Configuration)
    return(W_diag(e, c) + W_off_diag(e,c))
end

function E(e::Ensemble, c::Configuration)
    return(Ekin(e, c) + Epot(e,c))
end

function K(e::Ensemble, c::Configuration)
    return(length(c.kinks))
end

function occupations(e::Ensemble, c::Configuration, emax::Int=100)
    nk = zeros(emax)
    if isempty(c.kinks)
        for en in get_energy.(c.occupations)
            if en < emax
                nk[en+1] = nk[en+1] + 1
            end
        end
    else
        occs = copy(c.occupations)
        old_Tau = first(last(c.kinks)) - 1
        for (tau,k) in c.kinks
            for en in get_energy.(occs)
                if en < emax
                    nk[en+1] = nk[en+1] + (tau - old_Tau)
                end
            end
            old_Tau = tau
            change_occupations(occs,k)
        end
    end
    nk
end

function particleNumber(c::Configuration)
  return length(c.occupations)
end
