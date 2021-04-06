export Model

"""
    abstract type Model end

Abstract type for the representation of a physical model for which calculations can be performed.
Needs to be implemented system-specifically, in particular the following functions need to be defined: (for now)

- the one-particle matrix element, `energy(::Model, ::Orbital)`
- the two-particle matrix element, `w(::Model, ::Orbital, ::Orbital, ::Orbital, ::Orbital)`

"""
abstract type Model end

function energy(m::Model, o::Orbital)
    error("missing implementation of energy(::$(typeof(m)), ::$(typeof(o)))")
end

function w(m::Model, o1::Orbital, o2::Orbital, o3::Orbital, o4::Orbital)
    error("missing implementation of w(::$(typeof(m)), ::$(typeof(o1)), ::$(typeof(o2)))")
end

@doc raw"""
    wminus(::Model, i,j,k,l)

return the difference of the two-particle matrix element with the same but the last two indices transposed
antisymmetric difference of the two-particle matrix elements:

$w^-_{ijkl} = w_{ijkl} - w_{ijlk}$

This is called the antisymmetrized two-particle matrix element.
This is an abbreviation for these terms arising in the
Slater-Condon rules for the calculation of the many-body matrix elements
via the one- and two-particle matrix elements of the underlying single-particle basis.
"""
wminus(m::Model, i,j,k,l) = w(m,i,j,k,l) - w(m,i,j,l,k)


@doc raw"""
    Woffdiag_element(::Model, ::Ensemble, ::Orbital, ::Orbital, ::Orbital, ::Orbital)

Return the offdiagonal many body matrix element of the interaction operator

$\frac{1}{4} ( w_{ijkl} - w_{ijlk} ) (\pm \langle\{\tilde{n}\}|\{n\}^{ij}_{kl}\rangle)$

for an excitation given by creating orbitals i,j and annihilating orbitals k, l
as given by the Slater-Condon rules.
"""
function Woffdiag_element(m::Model, e::Ensemble, i::Orbital, j::Orbital, k::Orbital, l::Orbital)
    # We sample with the weight of antisymmetrized matrix element but we do not restrict
    # the order of indices of our possible kinks. We therefor need an extra factor 1/4 in the weight-function
    1/4 * e.λ * wminus(m, i,j,k,l)
end

Woffdiag_element(m::Model, e::Ensemble, kink::T4) = Woffdiag_element(m, e, kink.i, kink.j, kink.k, kink.l)

"""
    ΔWoffdiag_element(::Model, e::Ensemble, itr, itr)

return the change in the offdiagonal many body matrix element of the interaction operator
given by adding the kinks in the first iterable and removing the kinks in the second iterable

both arguments `itr` are required to be iterables containing kinks
"""
function ΔWoffdiag_element(m::Model, e::Ensemble, add_kinks, drop_kinks)
    if isempty(drop_kinks)
        prod(Woffdiag_element(m, e, k) for k in add_kinks)
    elseif isempty(add_kinks)
        1.0 / prod(Woffdiag_element(m, e, k) for k in drop_kinks)
    else
        prod(Woffdiag_element(m, e, k) for k in add_kinks) / prod(Woffdiag_element(m, e, k) for k in drop_kinks)
    end
end

ΔWoffdiag_element(m::Model, e::Ensemble, add_kinks::SortedDict{ImgTime,<:Kink}, drop_kinks::SortedDict{ImgTime,<:Kink}) = ΔWoffdiag_element(m, e, values(add_kinks), values(drop_kinks))


"""
    ΔW_diag(m::Model, i, j, k, l, occ)

change in the diagonal interaction matrix element due to a change in the occupation `occ`
in a periodic interval (τ1,τ2) where no kinks occur
the change in the occupation is assumed to consist in a creation of two orbitals i, j
and in the annihlation of two orbitals k, l
"""
function ΔW_diag(m::Model, i, j, k, l, occ)
    @assert (i ∉ occ) & (j ∉ occ) "Calculation of the change in the many-body diagonal interaction matrix element: This function assumes that the first two orbitals\n\t $(i)\n and\n\t $(j) given are the creator orbitals and thus that they are not occupied in the given occupation\n\t $(occ). "
    # contributions due to mean field interactions of the annihilated orbitals
    Δ = sum( wminus(m, ν,k,k,ν) + wminus(m, ν,l,l,ν) for ν in drop(occ, Set([k,l])) )# interactions of mean field with k and l
    Δ += wminus(m, k,l,l,k)# interaction between k and l
    # contributions due to mean field interactions of the created orbitals
    # note: the annihilator orbitals k, l are not in the new occupation
    Δ -= sum( wminus(m, ν,i,i,ν) + wminus(m, ν,j,j,ν) for ν in drop(occ, Set([k,l])) )# interactions of mean field with i and j
    Δ -= wminus(m, i,j,j,i)# interaction between i and j
    return Δ
end

"""
    ΔW_diag(m::Model, i, j, occ)

change in the diagonal interaction matrix element due to a change in the occupation `occ`
in a periodic interval (τ1,τ2) where no kinks occur
the change in the occupation is assumed to consist in a creation of one orbitals i
and in the annihlation of one orbitals j
"""
function ΔW_diag(m::Model, i, j, occ)
    @assert (i ∉ occ) "Calculation of the change in the many-body diagonal interaction matrix element: This function assumes that the first two orbitals\n\t $(i)\n and\n\t $(j) given are the creator orbitals and thus that they are not occupied in the given occupation\n\t $(occ). "
    # contributions due to mean field interactions of the annihilated orbitals
    Δ = sum( wminus(m, ν,j,j,ν) for ν in drop(occ, j) )# interactions of mean field with k
    # contributions due to mean field interactions of the created orbitals
    Δ -= sum( wminus(m, ν,i,i,ν) for ν in drop(occ, j) )
end



### convention: all ΔX_element represent only the difference in the matrix elements will be used as exp(-Δ) for the weight change
###             thus, contributions from the new (proposed) configuration (creators) appear positive (+)
###             and the contributions from the old configuration (annihilators) appear negative (-)

"""
    ΔT_element(::Model, i,j,k,l)

return the change in the kinetic many body matrix element due
to creating orbitals i, j and annihilating orbitals k, l
"""
ΔT_element(m::Model, i,j,k,l) = energy(m, i) + energy(m, j) - energy(m, k) - energy(m, l)



"""
    ΔWdiag_element(::Model, ::Ensemble, ::Configuration, i, j, k, l, τ1, τ2)

Calculate the change in the diagonal interaction many-body matrix element
due to a change in the occupations given by creating two orbitals i and j
and annihilating two orbitals k, l in the interval (τ1, τ2).
This interval may be periodically extended over the bounds (0,1) if τ1 > τ2,
i.e. the change in the occupation is considered for (τ1,1] ∪ [0,τ2) in that case.
We do not need to evaluate the diagonal interaction
between all orbitals in all time-intervalls, but it is sufficient to evaluate
the full diagonal interaction with the occupations at the start of the intervall
and then consider only contributions of orbitals that are changed by kinks in the intervall.
"""
function ΔWdiag_element(m::Model, e::Ensemble, c::Configuration, i, j, k, l, τ1, τ2)# TODO: assuming that i, j are creators and k, l are annihilators. Use Step instead ?
    @assert τ1 != τ2 " The diagonal interaction matrix element changes when kinks are added at different times and thus the occupations between the kinks are altered. It has no meaning to calculate this matrix element (or to add kinks) at equal times τ1=$(τ1), τ2=$(τ2). "

    # get all kinks between τ1 and τ2
    Ks = kinks_from_periodic_interval(c.kinks, τ1, τ2)

    # calculate Wdiag with the occupation at the start of the interval
    ΔWdiag = ΔW_diag(m, i,j,k,l, occupations(c,τ1)) * Δ(τ1,τ2)

    if !isempty(Ks)
        # calculate contrubutions to ΔWdiag from the orbitals changed by kinks in the intervall:
        # add a contribution if an orbital is created and
        # remove a contribution if an orbital is annilated
        # this is more efficient than calculating the occupation for each consecutive time-interval
        # via `occupation(occ,t)` since this function applies all kinks up to t::ImgTime
        τs = times_from_periodic_interval(c.kinks, τ1, τ2)# get a time-ordered list of the times of the kinks between τ1 and τ2
        ΔWdiag += sum( (ΔW_diag(m, i,j,k,l, creators(Ks[t1])) - ΔW_diag(m, i,j,k,l, annihilators(Ks[t1]))) * Δ(t1,τ2) for t1 in τs)
    end
    e.λ * ΔWdiag
end


"""
    ΔWdiag_element(::Model, ::Ensemble, ::Configuration, i, j, τ1, τ2)

Calculate the change in the diagonal interaction many-body matrix element
due to a change in the occupations given by creating an orbital i
and annihilating an orbital j in the interval (τ1, τ2).
This interval may be periodically extended over the bounds (0,1) if τ1 > τ2,
i.e. the change in the occupation is considered for (τ1,1] ∪ [0,τ2) in that case.
We do not need to evaluate the diagonal interaction
between all orbitals in all time-intervalls, but it is sufficient to evaluate
the full diagonal interaction with the occupations at the start of the intervall
and then consider only contributions of orbitals that are changed by kinks in the intervall.
"""
function ΔWdiag_element(m::Model, e::Ensemble, c::Configuration, i, j, τ1, τ2)# TODO: assuming that i is creator and j is annihilator. Use Step instead ?
    @assert τ1 != τ2 " The diagonal interaction matrix element changes when kinks are added at different times and thus the occupations between the kinks are altered. It has no meaning to calculate this matrix element (or to add kinks) at equal times τ1=$(τ1), τ2=$(τ2). "

    # get all kinks between τ1 and τ2
    Ks = kinks_from_periodic_interval(c.kinks, τ1, τ2)

    # calculate Wdiag with the occupation at the start of the interval
    ΔWdiag = ΔW_diag(m, i,j, occupations(c,τ1)) * Δ(τ1,τ2)
    if !isempty(Ks)
        # calculate contrubutions to ΔWdiag from the orbitals changed by kinks in the intervall:
        # add a contribution if an orbital is created and
        # remove a contribution if an orbital is annilated
        # this is more efficient than calculating the occupation for each consecutive time-interval
        # via `occupation(occ,t)` since this function applies all kinks up to t::ImgTime
        τs = times_from_periodic_interval(c.kinks, τ1, τ2)# get a time-ordered list of the times of the kinks between τ1 and τ2
        ΔWdiag += sum( (ΔW_diag(m, i,j, creators(Ks[t1])) - ΔW_diag(m, i,j, annihilators(Ks[t1]))) * Δ(t1,τ2) for t1 in τs)
    end
    e.λ * ΔWdiag
end


@doc raw"""
    sign_offdiagonal_product(::Model, ::Configuration)

Return the sign of the product of two-particle terms in the offdiagonal many-body matrix elements.
These are given by function `wminus`, i.e.

$\langle \{\tilde{n}\} | a^{\dagger}_i a^{\dagger}_j a_k a_l | \{n\} \rangle
    = \pm ( w_{ijkl} - w_{ijlk} ) \text{ for } \{\tilde{n}\} = \{n\}_{kl}^{ij}$

the term in the braces may be negative and this function returns the product of the
sign of all these contributions $w_{ijkl} - w_{ijlk}$ from all kinks.
The sign $\pm$ is determined from the permutation factor of the orbitals i,j,k,l
and is not calculated here.

used for the calculation the sign of the weight function
"""
function sign_offdiagonal_product(m::Model, c::Configuration)
    s = 1
    for κ in values(c.kinks)
        s *= sign(wminus(m, κ.i, κ.j, κ.k, κ.l))
    end
    return s
end

"""
    signum(m::Model, c::Configuration)

Calculate the sign of the `Configuration`'s weight.
"""
signum(m::Model, c::Configuration) = ladder_operator_order_factor(c.kinks)*sign_offdiagonal_product(m, c)
