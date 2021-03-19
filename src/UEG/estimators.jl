#include("../../../src/HEG/model.jl")

function Ekin(e::Ensemble, c::Configuration)
    if length(c.kinks) == 0
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

function W_off_diag(e::Ensemble, c::Configuration)
    return(-(length(c.kinks)/e.β))
end

function W_diag(e::Ensemble, c::Configuration)
    W_diag = 0
    if length(c.kinks) == 0
        for occ1 in c.occupations
            redundant = true
            for occ2 in c.occupations
                if (!redundant & (occ1.spin == occ2.spin))
                    W_diag += λ(e.N,e.rs) / dot((occ1.vec-occ2.vec),(occ1.vec-occ2.vec))
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
                        W_diag += λ(e.N,e.rs) / dot((occ1.vec-occ2.vec),(occ1.vec-occ2.vec)) * (τ-old_τ)
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

signum(e::Ensemble, c::Configuration) = signum(c)
signum(c::Configuration) = ladder_operator_order_factor(c.kinks)*sign_offdiagonal_product(c)

function particleNumber(c::Configuration)
  return length(c.occupations)
end

function W_off_diag(e::Ensemble, avg_K::Float64)
    return (-(avg_K/e.β))
end


#####################Calculationg observables after Simulation
function abs_E_mad(N::Int, lam::Float64) #internal units
    return 2.83729747948527 * pi/2.0 * N * lam
end

function E_int_from_Hartree(E_Ha::Float64, lam::Float64)
    return (E_Ha /(16/((2*pi)^4 * (lam/2)^2) * 0.5))
end

function E_int_from_Hartree(E_Ha::Float64, e::Ensemble)
    return (E_Ha /(16/((2*pi)^4 * (λ(e.N, e.rs)/2)^2) * 0.5))
end

function E_Ry(E_internal::Float64, lam::Float64)
    return (E_internal * 16/((2*pi)^4 * lam^2))
end

function E_Ry(E_internal::Float64, e::Ensemble)
    return (E_internal * 16/((2*pi)^4 * λ(e.N, e.rs)^2))
end

function E_Ha(E_internal::Float64, lam::Float64)
    return (E_internal * 16/((2*pi)^4 * lam^2) * 0.5)
end

function E_Ha(E_internal::Float64, e::Ensemble)
    return (E_internal * 16/((2*pi)^4 * λ(e.N, e.rs)^2) * 0.5)
end


function Et_Ry(E_internal::Float64, e::Ensemble)
    return (E_Ry(E_internal-abs_E_mad(e.N, λ(e.N,e.rs)),λ(e.N,e.rs)))
end
