export Model

"""
    abstract type Model end

Abstract type for the representation of a physical model for which calculations can be performed.
Needs to be implemented system-specifically, in particular the following functions need to be defined: (for now)

    * the one-particle matrix element, `energy(::Model, ::Orbital)`
    * interaction kernel, `kernel(::Model, ::Orbital, ::Orbital)`
"""
abstract type Model end

function energy(m::Model, o::Orbital)
    error("missing implementation of energy(::$(typeof(m)), ::$(typeof(o)))")
end

function kernel(m::Model, o1::Orbital, o2::Orbital)
    error("missing implementation of kernel(::$(typeof(m)), ::$(typeof(o1)), ::$(typeof(o2)))")
end


" anti-symmetrized interaction matrix element "
function wminus(m::Model, i::Orbital, j::Orbital, k::Orbital, l::Orbital) where {D}
    @assert ((i != k) && (i != l))
    @assert ((i.spin == j.spin) == (k.spin == l.spin))
    @assert(in(i.spin,[k.spin,l.spin]))
    @assert(i.vec + j.vec == k.vec + l.vec)
    if i.spin == j.spin
        @assert(!isinf(abs(kernel(m, i, k) - kernel(m, i, l))))
        return kernel(m, i, k) - kernel(m, i, l)
    elseif i.spin == k.spin
        @assert(!isinf(abs(kernel(m, i, k))))
        return kernel(m, i, k)
    elseif i.spin == l.spin
        @assert(!isinf(abs(kernel(m, i, l))))
        return -kernel(m, i, l)
    end
end

" diagonal interaction matrix element "
function wdiag(m::Model, a::Orbital, b::Orbital)
    if a.spin != b.spin
        return 0
    else
        return -kernel(m, a, b)
    end
end


wminus(m::Model, kink::T4) = wminus(m, kink.i, kink.j, kink.k, kink.l)




function offdiagonal_element(m::Model, e::Ensemble, kink::T4)
    # We sample with the weight of antisymmetrized matrix element but we do not restrict
    # the order of indices of our possible kinks. We therefor need an extra factor 1/4 in the weight-function
    return 1/4 * e.λ * wminus(m, kink)
end


"""
    Δdiagonal_interaction(::Model, ::Ensemble, ::Configuration, a::Orbital, b::Orbital, c::Orbital, d::Orbital, τ1, τ2)

change in the diagonal part of the interaction when changing the occupation by creating orbitals `a, b` and annihilating `c, d` between τ1 and τ2
"""
function Δdiagonal_interaction(m::Model, e::Ensemble, c::Configuration, orb_a::Orbital, orb_b::Orbital, orb_c::Orbital, orb_d::Orbital, τ1, τ2)
    Δτ12 = τ2 - τ1

    if Δτ12 < 0
        Δτ12 += 1
    end

    @assert (orb_a.spin != orb_b.spin) == (orb_c.spin != orb_d.spin )

    Δdi = Δτ12 * e.λ * ( wdiag(m, orb_a, orb_b) - wdiag(m, orb_c, orb_d) )

    occs = occupations(c, τ1)
    @assert !isinf(abs(Δdi))
    # collect diagonal interaction energy at τ1
    for occ in occs
        if occ in [ orb_a, orb_b, orb_c, orb_d ]
            continue
        else
            for orb in [orb_a, orb_b]
                Δdi += Δτ12 * e.λ * wdiag(m, occ, orb)
            end
            for orb in [orb_c, orb_d]
                Δdi -= Δτ12 * e.λ * wdiag(m, occ, orb)
            end
        end
        @assert !isinf(abs(Δdi))
        @assert(!isnan(Δdi))
    end

    if isempty(c.kinks)
        return Δdi
    end

    kink_semi_token = searchsortedfirst(c.kinks,τ1)
    if kink_semi_token == pastendsemitoken(c.kinks)
        kink_semi_token = startof(c.kinks)
    end
    τ_kink,kink = deref((c.kinks,kink_semi_token))

    # the kink at τ1 is already considered in occs
    if τ_kink == τ1
        kink_semi_token = advance((c.kinks,kink_semi_token))
        if kink_semi_token == pastendsemitoken(c.kinks)
            kink_semi_token = startof(c.kinks)
        end
        τ_kink,kink = deref((c.kinks,kink_semi_token))
    end

    loop_counter = 0

    @assert !isinf(abs(Δdi))
    @assert(!isnan(Δdi))
    # collect contributions to diagonal interaction energy due to kinks in the interval
    while ((τ1 < τ_kink < τ2) | (τ_kink < τ2 < τ1) | (τ2 < τ1 < τ_kink)) & (loop_counter < length(c.kinks))
        Δτ = τ2 - τ_kink
        if Δτ < 0
            Δτ += 1
        end
        for occ in [kink.i, kink.j]
            for orb in [orb_a, orb_b]
                 Δdi += Δτ * e.λ * wdiag(m, occ, orb)
            end
            for orb in [orb_c, orb_d]
                 Δdi -= Δτ * e.λ * wdiag(m, occ, orb)
            end
        end
        for occ in [kink.k, kink.l]
            for orb in [orb_a, orb_b]
                 Δdi -= Δτ * e.λ * wdiag(m, occ, orb)
            end
            for orb in [orb_c, orb_d]
                 Δdi += Δτ * e.λ * wdiag(m, occ, orb)
            end
        end

        @assert !isinf(abs(Δdi))
        @assert(!isnan(Δdi))
        kink_semi_token = advance((c.kinks,kink_semi_token))
        if kink_semi_token == pastendsemitoken(c.kinks)
            kink_semi_token = startof(c.kinks)
        end
        τ_kink,kink = deref((c.kinks,kink_semi_token))
        loop_counter += 1
    end
    @assert(!isnan(Δdi))
    return Δdi
end

" This function assumes for now that all kinks in c are of type 4. "
function Δdiagonal_interaction(m::Model, e::Ensemble, c::Configuration, orb_a::Orbital, orb_b::Orbital, τ1, τ2)
    Δτ12 = τ2 - τ1
    if Δτ12 < 0
        Δτ12 += 1
    end
    Δdi = 0
    occs = occupations(c, τ1)
    @assert !in(orb_a, occs)
    for occ in occs
        if occ.vec in [ orb_a.vec, orb_b.vec ]
            continue
        else
            Δdi += Δτ12 * e.λ * wdiag(m, occ, orb_a)
            Δdi -= Δτ12 * e.λ * wdiag(m, occ, orb_b)
        end
    end
    if isempty(c.kinks)
        return Δdi
    end
    kink_semi_token = searchsortedfirst(c.kinks,τ1)
    if kink_semi_token == pastendsemitoken(c.kinks)
        kink_semi_token = startof(c.kinks)
    end
    τ_kink,kink = deref((c.kinks,kink_semi_token))

    # The kink at τ1 is already considered in occs.
    if τ_kink == τ1
        kink_semi_token = advance((c.kinks,kink_semi_token))
        if kink_semi_token == pastendsemitoken(c.kinks)
            kink_semi_token = startof(c.kinks)
        end
        τ_kink,kink = deref((c.kinks,kink_semi_token))
    end

    loop_counter = 0

    while ((τ1 < τ_kink < τ2) | (τ_kink < τ2 < τ1) | (τ2 < τ1 < τ_kink)) & (loop_counter < length(c.kinks))
        Δτ = τ2 - τ_kink
        if Δτ < 0
            Δτ += 1
        end

        for occ in [kink.i, kink.j]
            Δdi += Δτ * e.λ * wdiag(m, occ, orb_a)
            Δdi -= Δτ * e.λ * wdiag(m, occ, orb_b)
        end
        for occ in [kink.k, kink.l]
            Δdi += Δτ * e.λ * wdiag(m, occ, orb_b)
            Δdi -= Δτ * e.λ * wdiag(m, occ, orb_a)
        end


        kink_semi_token = advance((c.kinks,kink_semi_token))
        if kink_semi_token == pastendsemitoken(c.kinks)
            kink_semi_token = startof(c.kinks)
        end
        τ_kink,kink = deref((c.kinks,kink_semi_token))
        loop_counter += 1
    end

    return Δdi
end

"""
    sign_offdiagonal_product(m::Model, c::Configuration)

sign of the product of the offdiagonal elements
used for the calculation the sign of the weight function
"""
function sign_offdiagonal_product(m::Model, c::Configuration)
    sign_ofd_prod = 1
    for kink in values(c.kinks)
        sign_ofd_prod *= sign(wminus(m, kink))
    end
    return sign_ofd_prod
end

"""
    signum(m::Model, c::Configuration)

Calculate the sign of the `Configuration`'s weight.
"""
signum(m::Model, c::Configuration) = ladder_operator_order_factor(c.kinks)*sign_offdiagonal_product(m, c)
