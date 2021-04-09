"""
generic estimators for the CPIMC formalism
"""
module Estimators

using ..CPIMC

import ..CPIMC: right_type_1_count, longest_type_1_chain_length, energy, w, signum

export E, Ekin, W, W_off_diag, W_diag, K, occupations, signum, longest_type_1_chain_length, right_type_1_count

"""
    Ekin(m::Model, e::Ensemble, c::Configuration)

estimator for the kinetic energy
    """
function Ekin(m::Model, e::Ensemble, c::Configuration)
    if isempty(c.kinks)
        return sum(energy(m, n) for n in c.occupations)
    end

    occs = copy(c.occupations)
    E_kin = 0
    # betrachte das komplette Intervall zwischen letztem und erstem Kink als
    # erstes intervall
    old_τ = first(last(c.kinks)) - 1
    for (τ,kink) in c.kinks
        E_kin += sum(energy(m, n) for n in occs) * float(τ - old_τ)
        old_τ = τ
        excite!(occs, kink)
    end
    return E_kin
end

"""
    W_off_diag(m::Model, e::Ensemble, c::Configuration)

estimator for the offdiagonal contribution to the interaction energy
This estimator is redundant because we can calculate the offdiagonal interaction energy from K_fermion.
"""
function W_off_diag(m::Model, e::Ensemble, c::Configuration)
    return -(length(c.kinks)/e.β)
end

"""
    W_diag(m::Model, e::Ensemble, c::Configuration)

estimator for the diagonal contribution to the interaction energy
"""
function W_diag(m::Model, e::Ensemble, c::Configuration)
    W_diag = 0

    if isempty(c.kinks)
        for occ1 in c.occupations
            redundant = true
            for occ2 in c.occupations
                if !redundant & (occ1.spin == occ2.spin)
                    W_diag += e.λ * w(m, occ1, occ2, occ2, occ1)
                end
                if occ1 == occ2
                    redundant = false
                end
            end
        end
    else
        occs = copy(c.occupations)
        old_τ = first(last(c.kinks)) - 1
        for (τ,kink) in c.kinks
            for occ1 in occs
                redundant = true
                for occ2 in occs
                    if !redundant & (occ1.spin == occ2.spin)
                        W_diag += e.λ * w(m, occ1, occ2, occ2, occ1) * (τ-old_τ)
                    end
                    if occ1 == occ2
                        redundant = false
                    end
                end
            end
            old_τ = τ
            excite!(occs, kink)
        end
    end
    W_diag *= -1
    return W_diag
end

"""
    W(m::Model, e::Ensemble, c::Configuration)

estimator for the interaction energy
This estimator is redundant because we can calculate the interaction energy from K_fermion and W_diag.
"""
function W(m::Model, e::Ensemble, c::Configuration)
    return W_diag(m, e, c) + W_off_diag(m, e,c)
end

"""
    E(m::Model, e::Ensemble, c::Configuration)

estimator for the interaction energy
This estimator is redundant because we can calculate the full energy as the sum of the kinetic and the interaction part.
"""
function E(m::Model, e::Ensemble, c::Configuration)
    return Ekin(m, e, c) + Epot(m, e, c)
end

"""
    K(m::Model, e::Ensemble, c::Configuration)

estimator for the number of kinks
"""
function K(m::Model, e::Ensemble, c::Configuration)
    return length(c.kinks)
end

"""
    occupations(m::Model, e::Ensemble, c::Configuration, emax::Int=100) :: Array{Float64,1}

estimator for the occupations of the emax:Int lowest single particle energy eigenvalues
"""
function occupations(m::Model, e::Ensemble, c::Configuration, emax::Int=100) :: Array{Float64,1}
    nk = zeros(Float64, emax)
    if isempty(c.kinks)
        for en in [ energy(m,n) for n in c.occupations ]
            if en < emax
                nk[en+1] = nk[en+1] + 1.0
            end
        end
    else
        occs = copy(c.occupations)
        old_τ = first(last(c.kinks)) - 1
        for (tau,k) in c.kinks
            for en in [ energy(m,n) for n in c.occupations ]
                if en < emax
                    nk[en+1] = nk[en+1] + float(tau - old_τ)
                end
            end
            old_τ = tau
            excite!(occs,k)
        end
    end
    nk
end

"""
    signum(m::Model, e::Ensemble, c::Configuration)

estimator for the sign of the weight function
"""
signum(m::Model, e::Ensemble, c::Configuration) = signum(m, c)


longest_type_1_chain_length(m::Model, e::Ensemble, c::Configuration) = longest_type_1_chain_length(c.kinks)

right_type_1_count(m::Model, e::Ensemble, c::Configuration) = right_type_1_count(c.kinks)

end
