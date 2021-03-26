"""
    measure(measurements, e, c)
calculate estimators and fit to measurements.
    """
function measure(measurements, e, c)
    s = signum(c)
    for (key,(stat,obs)) in measurements
        if in(key,[:sign, :K])
            fit!(stat, obs(e,c))
        else
            fit!(stat, obs(e,c)*s)
        end
    end
end

"""
    Ekin(e::Ensemble, c::Configuration)

estimator for the kinetic energy
    """
function Ekin(e::Ensemble, c::Configuration)
    if isempty(c.kinks)
        return(sum(energy(n) for n in c.occupations))
    end
    occs = copy(c.occupations)
    E_kin = 0
    # betrachte das komplette Intervall zwischen letztem und erstem Kink als
    # erstes intervall
    old_τ = first(last(c.kinks)) - 1
    for (τ,kink) in c.kinks
        E_kin += sum(energy(n) for n in occs) * float(τ - old_τ)
        old_τ = τ
        excite!(occs, kink)
    end
    return(E_kin)
end

"""
    W_off_diag(e::Ensemble, c::Configuration)

estimator for the offdiagonal contribution to the interaction energy
This estimator is redundant because we can calculate the offdiagonal interaction energy from K_fermion.
"""
function W_off_diag(e::Ensemble, c::Configuration)
    return -(length(c.kinks)/e.β)
end

"""
    W_diag(e::Ensemble, c::Configuration)

estimator for the diagonal contribution to the interaction energy
"""
function W_diag(e::Ensemble, c::Configuration)
    W_diag = 0
    if isempty(c.kinks)
        for occ1 in c.occupations
            redundant = true
            for occ2 in c.occupations
                if (!redundant & (occ1.spin == occ2.spin))
                    W_diag += λ(e.N,e.rs) * kernel(occ1,occ2)
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
                    if (!redundant & (occ1.spin == occ2.spin))
                        W_diag += λ(e.N,e.rs) * kernel(occ1,occ2) * (τ-old_τ)
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
    Epot(e::Ensemble, c::Configuration)

estimator for the interaction energy
This estimator is redundant because we can calculate the interaction energy from K_fermion and W_diag.
"""
function W(e::Ensemble, c::Configuration)
    return W_diag(e, c) + W_off_diag(e,c)
end

"""
    E(e::Ensemble, c::Configuration)

estimator for the interaction energy
This estimator is redundant because we can calculate the full energy as the sum of the kinetic and the interaction part.
"""
function E(e::Ensemble, c::Configuration)
    return Ekin(e, c) + Epot(e,c)
end

"""
    K(e::Ensemble, c::Configuration)

estimator for the number of kinks
"""
function K(e::Ensemble, c::Configuration)
    return length(c.kinks)
end

"""
    occupations(e::Ensemble, c::Configuration, emax::Int=100) :: Array{Float64,1}

estimator for the occupations of the emax:Int lowest single particle energy eigenvalues
"""
function occupations(e::Ensemble, c::Configuration, emax::Int=100) :: Array{Float64,1}
    nk = zeros(Float64, emax)
    if isempty(c.kinks)
        for en in energy.(c.occupations)
            if en < emax
                nk[en+1] = nk[en+1] + 1.0
            end
        end
    else
        occs = copy(c.occupations)
        old_τ = first(last(c.kinks)) - 1
        for (tau,k) in c.kinks
            for en in energy.(occs)
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
    signum(e::Ensemble, c::Configuration)

estimator for the sign of the weight function
"""
signum(e::Ensemble, c::Configuration) = signum(c)
signum(c::Configuration) = ladder_operator_order_factor(c.kinks)*sign_offdiagonal_product(c)


#####################Calculationg observables after Simulation
function W_off_diag(e::Ensemble, avg_K::Float64)
    return (-(avg_K/e.β))
end

function abs_E_madelung(N::Int, λ::Float64) #internal units
    return 2.83729747948527 * pi/2.0 * N * λ
end

function E_int_from_Hartree(E_Ha::Float64, λ::Float64)
    return (E_Ha /(16/((2*pi)^4 * (λ/2)^2) * 0.5))
end

function E_int_from_Hartree(E_Ha::Float64, e::Ensemble)
    return (E_Ha /(16/((2*pi)^4 * (λ(e.N, e.rs)/2)^2) * 0.5))
end

function E_Ry(E_internal::Float64, λ::Float64)
    return (E_internal * 16/((2*pi)^4 * λ^2))
end

function E_Ry(E_internal::Float64, e::Ensemble)
    return (E_internal * 16/((2*pi)^4 * λ(e.N, e.rs)^2))
end

function E_Ha(E_internal::Float64, λ::Float64)
    return (E_internal * 16/((2*pi)^4 * λ^2) * 0.5)
end

function E_Ha(E_internal::Float64, e::Ensemble)
    return (E_internal * 16/((2*pi)^4 * λ(e.N, e.rs)^2) * 0.5)
end


function Et_Ry(E_internal::Float64, e::Ensemble)
    return (E_Ry(E_internal-abs_E_madelung(e.N, λ(e.N,e.rs)),λ(e.N,e.rs)))
end

function Et_Ha(E_internal::Float64, e::Ensemble)
    return E_Ha(E_internal-abs_E_madelung(e.N, λ(e.N,e.rs)),λ(e.N,e.rs))
end
