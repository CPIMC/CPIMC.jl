
S = sphere_with_same_spin(OrbitalHEG((0,0,0)),dk=1)
a = OrbitalHEG((-2,0,0))
b = OrbitalHEG((3,0,0))
c = OrbitalHEG((0,0,0))
d = OrbitalHEG((1,0,0))
e = OrbitalHEG((5,9,9))
f = OrbitalHEG((1,1,1))

g = OrbitalHEG(a.vec + e.vec - f.vec, Up)
h = OrbitalHEG(c.vec + d.vec - g.vec, Up)



sd = SortedDict{ImgTime, Kink{<:Orbital}}( ImgTime(0.2) => T4(a,b,c,d),
                                           ImgTime(0.5) => T4(c,d,a,b),
                                           ImgTime(0.6) => T4(b,a,d,c),
                                           ImgTime(0.8) => T4(d,c,b,a) )

conf = Configuration(sphere(OrbitalHEG((0,0,0),Up),dk=1),sd)

@testset "orbs" begin
    @test orbs(T4(a,b,c,d)) == Set([a,b,c,d])
end

@testset "adjacent_kinks_affecting_orbs" begin
    @test adjacent_kinks_affecting_orbs(conf, Set([a]), ImgTime(0.0)) == (ImgTime(0.8) => T4(d,c,b,a), ImgTime(0.2) => T4(a,b,c,d))
    @test adjacent_kinks_affecting_orbs(conf, Set([b]), ImgTime(0.1)) == (ImgTime(0.8) => T4(d,c,b,a), ImgTime(0.2) => T4(a,b,c,d))
    @test adjacent_kinks_affecting_orbs(conf, Set([c]), ImgTime(0.2)) == (ImgTime(0.8) => T4(d,c,b,a), ImgTime(0.5) => T4(c,d,a,b))
    @test adjacent_kinks_affecting_orbs(conf, Set([d]), ImgTime(0.3)) == (ImgTime(0.2) => T4(a,b,c,d), ImgTime(0.5) => T4(c,d,a,b))
    @test adjacent_kinks_affecting_orbs(conf, Set([a,b]), ImgTime(0.4)) == (ImgTime(0.2) => T4(a,b,c,d), ImgTime(0.5) => T4(c,d,a,b))
    @test adjacent_kinks_affecting_orbs(conf, Set([b,e]), ImgTime(0.5)) == (ImgTime(0.2) => T4(a,b,c,d), ImgTime(0.6) => T4(b,a,d,c))
    @test adjacent_kinks_affecting_orbs(conf, Set([c,d]), ImgTime(0.6)) == (ImgTime(0.5) => T4(c,d,a,b), ImgTime(0.8) => T4(d,c,b,a))
    @test adjacent_kinks_affecting_orbs(conf, Set([d,a]), ImgTime(0.7)) == (ImgTime(0.6) => T4(b,a,d,c), ImgTime(0.8) => T4(d,c,b,a))
    @test adjacent_kinks_affecting_orbs(conf, Set([a,c,d]), ImgTime(0.8)) == (ImgTime(0.6) => T4(b,a,d,c), ImgTime(0.2) => T4(a,b,c,d))
    @test adjacent_kinks_affecting_orbs(conf, Set([b,c,d]), ImgTime(0.9)) == (ImgTime(0.8) => T4(d,c,b,a), ImgTime(0.2) => T4(a,b,c,d))
    @test adjacent_kinks_affecting_orbs(conf, Set([e]), ImgTime(0.9)) == (nothing,nothing)
end

@testset "kinks_affecting_orbs" begin
    @test kinks_affecting_orbs(conf, Set([a])) == sd
    @test kinks_affecting_orbs(conf, Set([b])) == sd
    @test kinks_affecting_orbs(conf, Set([c])) == sd
    @test kinks_affecting_orbs(conf, Set([d])) == sd
    @test kinks_affecting_orbs(conf, Set([e])) == SortedDict{ImgTime, Kink{<:Orbital}}()
end

@testset "τ_borders" begin
    @test τ_borders(conf, Set([a]), ImgTime(0.0)) == (ImgTime(0.8), ImgTime(0.2))
    @test τ_borders(conf, Set([a,e]), ImgTime(0.5)) == (ImgTime(0.2), ImgTime(0.6))
    @test τ_borders(conf, Set([b,d]), ImgTime(0.1)) == (ImgTime(0.8), ImgTime(0.2))
    @test τ_borders(conf, Set([e]), ImgTime(0.9)) == (ImgTime(0.0), ImgTime(1.0))
end

@testset "isunaffected" begin
    @test !isunaffected(conf.kinks, a)
    @test !isunaffected(conf.kinks, b)
    @test !isunaffected(conf.kinks, c)
    @test !isunaffected(conf.kinks, d)
    @test isunaffected(conf.kinks, e)
end

@testset "time_ordered_orbs(::T4)" begin
    @test time_ordered_orbs(T4(a,b,c,d))[1] == a
    @test time_ordered_orbs(T4(a,b,c,d))[2] == b
    @test time_ordered_orbs(T4(a,b,c,d))[3] == c
    @test time_ordered_orbs(T4(a,b,c,d))[4] == d
    @test time_ordered_orbs(T4(a,b,c,d)) == [a,b,c,d]
end

@testset "time_ordered_orbs(::SortedDict{ImgTime,<:Kink})" begin
    @test time_ordered_orbs(conf.kinks)[1] == a
    @test time_ordered_orbs(conf.kinks)[4] == d
    @test time_ordered_orbs(conf.kinks)[end] == a
    @test time_ordered_orbs(conf.kinks)[end-4] == c
end

@testset "kinks_from_periodic_interval(ck::SortedDict{ImgTime,<:Kink}, τ1, τ2)" begin
    ImgTime(0.2) => T4(a,b,c,d), ImgTime(0.5) => T4(c,d,a,b), ImgTime(0.6) => T4(b,a,d,c), ImgTime(0.8) => T4(d,c,b,a)
    @test kinks_from_periodic_interval(sd, 0.1, 0.6) == SortedDict(ImgTime(0.2) => T4(a,b,c,d), ImgTime(0.5) => T4(c,d,a,b))
    @test kinks_from_periodic_interval(sd, 0.7, 0.4) == SortedDict(ImgTime(0.2) => T4(a,b,c,d), ImgTime(0.8) => T4(d,c,b,a))
    @test kinks_from_periodic_interval(sd, 0.5, 0.7) == SortedDict(ImgTime(0.6) => T4(b,a,d,c))

    @test isempty(kinks_from_periodic_interval(sd, 0.3, 0.4))
    @test kinks_from_periodic_interval(sd, 0.4, 0.3) == sd

    @test isempty(kinks_from_periodic_interval(sd, 0.5, 0.5))
    @test isempty(kinks_from_periodic_interval(SortedDict{ImgTime,T4}(), 0.1, 0.9))
end

@testset "times_from_periodic_interval" begin
    t1 = ImgTime(0.1)
    t2 = ImgTime(0.3)
    @test times_from_periodic_interval(sd, t1, t2) == [t1, ImgTime(0.2), t2]
    @test times_from_periodic_interval(sd, t2, t1) == [t2,  ImgTime(0.5), ImgTime(0.6), ImgTime(0.8), t1]

    # test for times with no kinks in between
    t1 = ImgTime(0.3)
    t2 = ImgTime(0.4)
    @test times_from_periodic_interval(sd, t1, t2) == [t1,t2]
end

@testset "Δ(τ1::ImgTime,τ2::ImgTime)" begin
    @test Δ(ImgTime(0.5),ImgTime(0.3)) == ImgTime(0.5) - ImgTime(0.3)
    @test Δ(ImgTime(0.2),ImgTime(0.8)) == ImgTime(1) + ImgTime(0.2) - ImgTime(0.8)
    @test iszero( Δ(ImgTime(0.8),ImgTime(0.8)) )

    @test float(Δ(ImgTime(0.5),ImgTime(0.3))) ≈ float(ImgTime(0.2))
end

@testset "Type_1_investigation" begin
    g = OrbitalHEG(a.vec + e.vec - f.vec, Up)
    h = OrbitalHEG(c.vec + d.vec - g.vec, Up)

    Type_1_chain = SortedDict{ImgTime, Kink{<:Orbital}}( ImgTime(0.2) => T4(a,b,c,d),
                                               ImgTime(0.5) => T4(f,g,e,a),
                                               ImgTime(0.6) => T4(c,d,h,g),
                                               ImgTime(0.8) => T4(e,h,b,f) )

    @test (a.vec + b.vec - c.vec - d.vec) == OrbitalHEG((0,0,0)).vec
    @test (f.vec + g.vec - e.vec - a.vec) == OrbitalHEG((0,0,0)).vec
    @test (c.vec + d.vec - h.vec - g.vec) == OrbitalHEG((0,0,0)).vec
    @test (e.vec + h.vec - b.vec - f.vec) == OrbitalHEG((0,0,0)).vec

    occs = setdiff!(union!(sphere(OrbitalHEG((0,0,0),Up),dk=1),
                        Set([e,h])),
                Set([g,f]))
    conf_Type_1 = Configuration(occs,Type_1_chain)
    @test (occupations(conf_Type_1, ImgTime(0.9)) == occs)
    @test longest_type_1_chain_length(conf_Type_1.kinks) == 4
    @test right_type_1_count(conf_Type_1.kinks) == 4

    Type_1_chain[0.52] = T4(e,a,g,f)
    Type_1_chain[0.54] = T4(g,f,a,e)

    conf_Type_1 = Configuration(occs,Type_1_chain)

    @test length(conf_Type_1.kinks) == 6

    @test longest_type_1_chain_length(conf_Type_1.kinks) == 4
    @test right_type_1_count(conf_Type_1.kinks) == 4

    Type_1_chain[0.82] = T4(b,a,d,c)
    Type_1_chain[0.84] = T4(c,d,a,b)

    conf_Type_1 = Configuration(occs,Type_1_chain)
    @test longest_type_1_chain_length(conf_Type_1.kinks) == 3
    @test right_type_1_count(conf_Type_1.kinks) == 4

end


@testset "ΔWoffdiag_element" begin

    ens = CEnsemble(2.0, 5.680898543560106, 7)# θ: 0.125, λ: 0.09945178864947428

    t1 = ImgTime(0.3)
    t2 = ImgTime(0.7)
    τ3 = ImgTime(0.9)

    @test ΔWoffdiag_element(ens, kinks_from_periodic_interval(sd, ImgTime(0.3), ImgTime(0.9)), kinks_from_periodic_interval(sd, ImgTime(0.3), ImgTime(0.7))) ≈ Woffdiag_element(ens, d,c,b,a)
end






" coulomb kernel for 3D plane wavevectors "# for comparison of Δdiagonal_interaction
kernel(i::StaticVector{N,Int}, k::StaticVector{N,Int}) where {N} = 1.0 / dot( i-k, i-k )

" coulomb kernel in plane wave basis "# for comparison of Δdiagonal_interaction
kernel(i::OrbitalHEG,k::OrbitalHEG) = kernel(i.vec,k.vec)


" anti-symmetrized interaction matrix element "# for comparison of Δdiagonal_interaction
function wminus(i::Orbital, j::Orbital, k::Orbital, l::Orbital) where {D}
    @assert ((i != k) && (i != l))
    @assert ((i.spin == j.spin) == (k.spin == l.spin))
    @assert(in(i.spin,[k.spin,l.spin]))
    @assert(i.vec + j.vec == k.vec + l.vec)
    if i.spin == j.spin
        @assert(!isinf(abs(kernel(i, k) - kernel(i, l))))
        return kernel(i, k) - kernel(i, l)
    elseif i.spin == k.spin
        @assert(!isinf(abs(kernel(i, k))))
        return kernel(i, k)
    elseif i.spin == l.spin
        @assert(!isinf(abs(kernel(i, l))))
        return -kernel(i, l)
    end
end

" diagonal interaction matrix element "# for comparison Δdiagonal_interaction
function wdiag(a::Orbital, b::Orbital) where {D}
    if a.spin != b.spin
        return 0
    else
        return -kernel(a, b)
    end
end


wminus(kink::T4) = wminus(kink.i, kink.j, kink.k, kink.l)# for comparison of Δdiagonal_interaction


"""
    Δdiagonal_interaction(c::Configuration, e::Ensemble, orb_a::Orbital, orb_b::Orbital, orb_c::Orbital, orb_d::Orbital, τ1, τ2)

old implementation of this function
"""
function Δdiagonal_interaction(c::Configuration, e::Ensemble, orb_a::Orbital, orb_b::Orbital, orb_c::Orbital, orb_d::Orbital, τ1, τ2)
    Δτ12 = τ2 - τ1

    if Δτ12 < 0
        Δτ12 += 1
    end

    @assert (orb_a.spin != orb_b.spin) == (orb_c.spin != orb_d.spin )

    Δdi = Δτ12 * e.λ * ( wdiag(orb_a, orb_b) - wdiag(orb_c, orb_d) )

    occs = occupations(c, τ1)

    # collect diagonal interaction energy at τ1
    for occ in occs
        if occ in [ orb_a, orb_b, orb_c, orb_d ]
            continue
        else
            for orb in [orb_a, orb_b]
                Δdi += Δτ12 * e.λ * wdiag(occ,orb)
            end
            for orb in [orb_c, orb_d]
                Δdi -= Δτ12 * e.λ * wdiag(occ,orb)
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
                 Δdi += Δτ * e.λ * wdiag(occ,orb)
            end
            for orb in [orb_c, orb_d]
                 Δdi -= Δτ * e.λ * wdiag(occ,orb)
            end
        end
        for occ in [kink.k, kink.l]
            for orb in [orb_a, orb_b]
                 Δdi -= Δτ * e.λ * wdiag(occ,orb)
            end
            for orb in [orb_c, orb_d]
                 Δdi += Δτ * e.λ * wdiag(occ,orb)
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

" old implementation of this function "
function Δdiagonal_interaction(c::Configuration, e::Ensemble, orb_a::Orbital, orb_b::Orbital, τ1, τ2)
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
            Δdi += Δτ12 * e.λ * wdiag(occ,orb_a)
            Δdi -= Δτ12 * e.λ * wdiag(occ,orb_b)
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
            Δdi += Δτ * e.λ * wdiag(occ,orb_a)
            Δdi -= Δτ * e.λ * wdiag(occ,orb_b)
        end
        for occ in [kink.k, kink.l]
            Δdi += Δτ * e.λ * wdiag(occ,orb_b)
            Δdi -= Δτ * e.λ * wdiag(occ,orb_a)
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


@testset "Δdiagonal_interaction: compare with previous implementation" begin

    ens = CEnsemble(2.0, 5.680898543560106, 7)# θ: 0.125, λ: 0.09945178864947428

    # choose creator orbitals
    i = OrbitalHEG((1,1,1))
    j = OrbitalHEG((-1,-1,-1))
    # choose annihilator orbitals
    k = OrbitalHEG((0,-1,0))
    l = OrbitalHEG((0,1,0))
    @assert (i.spin == k.spin) & (j.spin == l.spin) " spin is not conserved for this excitation "
    @assert iszero( i.vec + j.vec - k.vec - l.vec ) " momentum is not conserved for this excitation "

    # chose times with no kink in between
    τ1 = ImgTime(0.3)
    τ2 = ImgTime(0.4)
    @test ΔWdiag_element(conf, ens, i, k, τ1, τ2) ≈ Δdiagonal_interaction(conf, ens, i, k, τ1, τ2)

    @assert all( Set([k,l]) .∈ (occupations(conf,τ1),) ) & all( Set([k,l]) .∈ (occupations(conf,τ2),) ) " the orbitals \n i=$i,\n j=$j,\n k=$k,\n l=$l cannot form a type B kink pair at times ($(float(τ1)), $(float(τ2))) "

    @test ΔWdiag_element(conf, ens, i, j, k, l, τ1, τ2) ≈ Δdiagonal_interaction(conf, ens, i, j, k, l, τ1, τ2)

    @test ΔWdiag_element(conf, ens, i, j, k, l, τ2, τ1) ≈ Δdiagonal_interaction(conf, ens, i, j, k, l, τ2, τ1)

    # choose times with one kink in between
    τ1 = ImgTime(0.1)
    τ2 = ImgTime(0.3)

    @assert all( Set([k,l]) .∈ (occupations(conf,τ1),) ) & all( Set([k,l]) .∈ (occupations(conf,τ2),) ) " the orbitals \n i=$i,\n j=$j,\n k=$k,\n l=$l cannot form a type B kink pair at times ($(float(τ1)), $(float(τ2))) "

    @test ΔWdiag_element(conf, ens, i, j, k, l, τ1, τ2) ≈ Δdiagonal_interaction(conf, ens, i, j, k, l, τ1, τ2)

    @test ΔWdiag_element(conf, ens, i, j, k, l, τ2, τ1) ≈ Δdiagonal_interaction(conf, ens, i, j, k, l, τ2, τ1)

    @test ΔWdiag_element(conf, ens, i, k, τ1, τ2) ≈ Δdiagonal_interaction(conf, ens, i, k, τ1, τ2)
    # Test for different spin case
    i = OrbitalHEG((1,1,1),Down)
    k = OrbitalHEG((0,-1,0), Down)
    @assert (i.spin == k.spin) & (j.spin == l.spin) " spin is not conserved for this excitation "
    @assert iszero( i.vec + j.vec - k.vec - l.vec ) " momentum is not conserved for this excitation "

    @test ΔWdiag_element(conf, ens, i, k, τ1, τ2) ≈ Δdiagonal_interaction(conf, ens, i, k, τ1, τ2)

    @test ΔWdiag_element(conf, ens, j, l, τ1, τ2) ≈ Δdiagonal_interaction(conf, ens, j, l, τ1, τ2)

    @test wdiag(i,j) == 0
    @test wdiag(k,l) == 0
    @test w(i,j,j,i) == 0
    @test w(k,l,l,k) == 0

    @test ΔWdiag_element(conf, ens, i, j, k, l, τ1, τ2) ≈ (ΔWdiag_element(conf, ens, i, k, τ1, τ2) + ΔWdiag_element(conf, ens, j, l, τ1, τ2))
    @test Δdiagonal_interaction(conf, ens, i, j, k, l, τ1, τ2) ≈ (ΔWdiag_element(conf, ens, i, k, τ1, τ2) + ΔWdiag_element(conf, ens, j, l, τ1, τ2))
    @test ΔWdiag_element(conf, ens, i, j, k, l, τ1, τ2) ≈ Δdiagonal_interaction(conf, ens, i, j, k, l, τ1, τ2)


    ### Test 1-particle excitation
    i = OrbitalHEG((-3,0,0))
    j = OrbitalHEG((0,1,0))

    @test ΔWdiag_element(conf, ens, i, j, τ1, τ2) ≈ Δdiagonal_interaction(conf, ens, i, j, τ1, τ2)

    @test ΔWdiag_element(conf, ens, i, j, τ2, τ1) ≈ Δdiagonal_interaction(conf, ens, i, j, τ2, τ1)


end


@testset "ΔWdiag_element(c::Configuration, e::Ensemble, i, j, k, l, τ1, τ2)" begin

    i = OrbitalHEG((0,-4,0))
    j = OrbitalHEG((0,3,1))
    k = OrbitalHEG((0,0,1))
    l = OrbitalHEG((0,-1,0))

    τ1 = ImgTime(0.1)
    τ2 = ImgTime(0.3)
    τ3 = ImgTime(0.4)

    λ = 0.8
    β = 0.02
    N = length(conf.occupations)

    @test ΔWdiag_element(conf, CEnsemble(λ, β, N), i, j, k, l, τ2, τ3) ≈ ΔW_diag(i, j, k, l, occupations(conf,τ2)) * (τ3 - τ2) * λ

    @test ΔWdiag_element(conf, CEnsemble(λ, β, N), i, j, k, l, τ1, τ3) ≈ ( ΔW_diag(i, j, k, l, occupations(conf,τ1)) * (ImgTime(0.2) - τ1)
                                                                                + ΔW_diag(i, j, k, l, occupations(conf,ImgTime(0.2))) * (τ3 - ImgTime(0.2))
                                                                                ) * λ

    @test ΔWdiag_element(conf, CEnsemble(λ, β, N), i, j, k, l, τ3, τ1) ≈ ( ΔW_diag(i, j, k, l, occupations(conf,τ3)) * (ImgTime(0.5) - τ3)
                                                                                + ΔW_diag(i, j, k, l, occupations(conf,ImgTime(0.5))) * (ImgTime(0.6) - ImgTime(0.5))
                                                                                + ΔW_diag(i, j, k, l, occupations(conf,ImgTime(0.6))) * (ImgTime(0.8) - ImgTime(0.6))
                                                                                + ΔW_diag(i, j, k, l, occupations(conf,ImgTime(0.8))) * (ImgTime(1) + τ1 - ImgTime(0.8))
                                                                                ) * λ
end
