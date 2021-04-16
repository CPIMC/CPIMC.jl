using CPIMC, CPIMC.PlaneWaves, CPIMC.UniformElectronGas, CPIMC.DefaultUpdates, DataStructures
import CPIMC: orbs, adjacent_kinks_affecting_orbs, kinks_affecting_orbs, prev, next, prev_affecting, next_affecting, τ_prev_affecting, τ_next_affecting, τ_borders, isunaffected, time_ordered_orbs, occupations_at, longest_type_1_chain_length, right_type_1_count, kinks_from_periodic_interval, times_from_periodic_interval, Δ, Woffdiag_element, ΔWoffdiag_element, ΔWdiag_element, ΔW_diag, add_orbs, add_orbs!, drop_orbs, drop_orbs!, drop_kinks, drop_kinks!, add_kinks, add_kinks!


S = sphere_with_same_spin(PlaneWave((0,0,0)),dk=1)
a = PlaneWave((-2,0,0))
b = PlaneWave((3,0,0))
c = PlaneWave((0,0,0))
d = PlaneWave((1,0,0))
e = PlaneWave((5,9,9))
f = PlaneWave((1,1,1))

g = PlaneWave(a.vec + e.vec - f.vec, Up)
h = PlaneWave(c.vec + d.vec - g.vec, Up)



sd = Kinks( ImgTime(0.2) => T4(a,b,c,d),
       ImgTime(0.5) => T4(c,d,a,b),
       ImgTime(0.6) => T4(b,a,d,c),
       ImgTime(0.8) => T4(d,c,b,a) )

conf = Configuration(sphere(PlaneWave((0,0,0),Up),dk=1),sd)

@testset "orbs" begin
    @test orbs(T4(a,b,c,d))[1] == a
    @test orbs(T4(a,b,c,d))[2] == b
    @test orbs(T4(a,b,c,d))[3] == c
    @test orbs(T4(a,b,c,d))[4] == d
    @test orbs(T4(a,b,c,d)) == (a,b,c,d)
    @test orbs(T2(a,b)) == (a,b)
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
    @test kinks_affecting_orbs(conf, Set([e])) == Kinks{PlaneWave}()
end

@testset "prev" begin
    @test prev(conf.kinks, ImgTime(0.0)) == 4
    @test prev(conf.kinks, ImgTime(0.5)) == 1
    @test prev(conf.kinks, ImgTime(0.1)) == 4
    @test prev(conf.kinks, ImgTime(0.7)) == 3
end

@testset "next" begin
    @test next(conf.kinks, ImgTime(0.1)) == 1
    @test next(conf.kinks, ImgTime(0.5)) == 3
    @test next(conf.kinks, ImgTime(0.0)) == 1
    @test next(conf.kinks, ImgTime(0.8)) == 1
end

@testset "prev_affecting" begin
    @test prev_affecting(conf.kinks, (a,), ImgTime(0.0)) == 4
    @test prev_affecting(conf.kinks, (a,e), ImgTime(0.5)) == 1
    @test prev_affecting(conf.kinks, (b,d), ImgTime(0.1)) == 4
    @test prev_affecting(conf.kinks, (b,d), ImgTime(0.1)) == 4
    @test prev_affecting(conf.kinks, (e,f), ImgTime(0.1)) == 0
end

@testset "next_affecting" begin
    @test next_affecting(conf.kinks, (a,), ImgTime(0.0)) == 1
    @test next_affecting(conf.kinks, (a,e), ImgTime(0.5)) == 3
    @test next_affecting(conf.kinks, (b,d), ImgTime(0.1)) == 1
    @test next_affecting(conf.kinks, (b,d), ImgTime(0.8)) == 1
    @test next_affecting(conf.kinks, (e,f), ImgTime(0.1)) == 0
end

@testset "τ_prev_affecting" begin
    @test τ_prev_affecting(conf.kinks, Set([a]), ImgTime(0.0)) == ImgTime(0.8)
    @test τ_prev_affecting(conf.kinks, Set([a,e]), ImgTime(0.5)) == ImgTime(0.2)
    @test τ_prev_affecting(conf.kinks, Set([b,d]), ImgTime(0.1)) == ImgTime(0.8)
end

@testset "τ_next_affecting" begin
    @test τ_next_affecting(conf.kinks, Set([a]), ImgTime(0.0)) == ImgTime(0.2)
    @test τ_next_affecting(conf.kinks, Set([a,e]), ImgTime(0.5)) == ImgTime(0.6)
    @test τ_next_affecting(conf.kinks, Set([b,d]), ImgTime(0.1)) == ImgTime(0.2)
end


@testset "isunaffected" begin
    @test !isunaffected(conf.kinks, a)
    @test !isunaffected(conf.kinks, b)
    @test !isunaffected(conf.kinks, c)
    @test !isunaffected(conf.kinks, d)
    @test isunaffected(conf.kinks, e)
end

@testset "time_ordered_orbs(::SortedDict{ImgTime,<:Kink})" begin
    @test time_ordered_orbs(conf.kinks)[1] == a
    @test time_ordered_orbs(conf.kinks)[4] == d
    @test time_ordered_orbs(conf.kinks)[end] == a
    @test time_ordered_orbs(conf.kinks)[end-4] == c
end

@testset "kinks_from_periodic_interval(ck::SortedDict{ImgTime,<:Kink}, τ1, τ2)" begin
    ImgTime(0.2) => T4(a,b,c,d), ImgTime(0.5) => T4(c,d,a,b), ImgTime(0.6) => T4(b,a,d,c), ImgTime(0.8) => T4(d,c,b,a)
    @test kinks_from_periodic_interval(sd, 0.1, 0.6) == Kinks(ImgTime(0.2) => T4(a,b,c,d), ImgTime(0.5) => T4(c,d,a,b))
    @test kinks_from_periodic_interval(sd, 0.7, 0.4) == Kinks(ImgTime(0.2) => T4(a,b,c,d), ImgTime(0.8) => T4(d,c,b,a))
    @test kinks_from_periodic_interval(sd, 0.5, 0.7) == Kinks(ImgTime(0.6) => T4(b,a,d,c))

    @test isempty(kinks_from_periodic_interval(sd, 0.3, 0.4))
    @test kinks_from_periodic_interval(sd, 0.4, 0.3) == sd

    @test isempty(kinks_from_periodic_interval(sd, 0.5, 0.5))
    @test isempty(kinks_from_periodic_interval(Kinks{PlaneWave{3}}(), 0.1, 0.9))
end

@testset "times_from_periodic_interval" begin
    t1 = ImgTime(0.1)
    t2 = ImgTime(0.3)
    @test times_from_periodic_interval(sd, t1, t2) == [ImgTime(0.2)]
    @test times_from_periodic_interval(sd, t2, t1) == [ImgTime(0.5), ImgTime(0.6), ImgTime(0.8)]

    # test for times with no kinks in between
    t1 = ImgTime(0.3)
    t2 = ImgTime(0.4)
    @test times_from_periodic_interval(sd, t1, t2) == []
end

@testset "Δ(τ1::ImgTime,τ2::ImgTime)" begin
    @test Δ(ImgTime(0.3), ImgTime(0.5)) == ImgTime(0.5) - ImgTime(0.3)
    @test Δ(ImgTime(0.8), ImgTime(0.2)) == ImgTime(1) + ImgTime(0.2) - ImgTime(0.8)
    @test iszero( Δ(ImgTime(0.8),ImgTime(0.8)) )

    @test float(Δ(ImgTime(0.3), ImgTime(0.5))) ≈ float(ImgTime(0.2))
end

@testset "Type_1_investigation" begin
    g = PlaneWave(a.vec + e.vec - f.vec, Up)
    h = PlaneWave(c.vec + d.vec - g.vec, Up)

    Type_1_chain = Kinks( ImgTime(0.2) => T4(a,b,c,d),
                          ImgTime(0.5) => T4(f,g,e,a),
                          ImgTime(0.6) => T4(c,d,h,g),
                          ImgTime(0.8) => T4(e,h,b,f) )

    @test (a.vec + b.vec - c.vec - d.vec) == PlaneWave((0,0,0)).vec
    @test (f.vec + g.vec - e.vec - a.vec) == PlaneWave((0,0,0)).vec
    @test (c.vec + d.vec - h.vec - g.vec) == PlaneWave((0,0,0)).vec
    @test (e.vec + h.vec - b.vec - f.vec) == PlaneWave((0,0,0)).vec

    occs = setdiff!(union!(sphere(PlaneWave((0,0,0),Up),dk=1), Set([e,h])), Set([g,f]))
    conf_Type_1 = Configuration(occs,Type_1_chain)
    @test (occupations_at(conf_Type_1, ImgTime(0.9)) == occs)
    @test longest_type_1_chain_length(conf_Type_1.kinks) == 4
    @test right_type_1_count(conf_Type_1.kinks) == 4

    add_kinks!(Type_1_chain, (ImgTime(0.52) => T4(e,a,g,f), ImgTime(0.54) => T4(g,f,a,e)))

    conf_Type_1 = Configuration(occs,Type_1_chain)

    @test length(conf_Type_1.kinks) == 6

    @test longest_type_1_chain_length(conf_Type_1.kinks) == 4
    @test right_type_1_count(conf_Type_1.kinks) == 4

    add_kinks!(Type_1_chain, (ImgTime(0.82) => T4(b,a,d,c), ImgTime(0.84) => T4(c,d,a,b)))

    conf_Type_1 = Configuration(occs,Type_1_chain)
    @test longest_type_1_chain_length(conf_Type_1.kinks) == 3
    @test right_type_1_count(conf_Type_1.kinks) == 4

end

@testset "find_fourth_orb_for_kink" for _ in (1:5)
    r = (1:100)
    orb1 = PlaneWave((rand(r),rand(r),rand(r)),rand([Up,Down]))
    orb2 = orb1
    while orb1 == orb2
        orb2 = PlaneWave((rand(r),rand(r),rand(r)),rand([Up,Down]))
    end
    orb3 = orb1
    while (orb3 == orb2) || (orb3 == orb1)
        orb3 = PlaneWave((rand(r),rand(r),rand(r)),rand([orb1.spin,orb2.spin]))
    end
    @test in(find_fourth_orb_for_kink(orb3, orb1, orb2).spin, [orb1.spin, orb2.spin])
    @test (find_fourth_orb_for_kink(orb3, orb1, orb2).spin == orb3.spin) == (orb1.spin == orb2.spin)
    @test find_fourth_orb_for_kink(orb3, orb1, orb2).vec + orb3.vec == orb1.vec + orb2.vec
end

@testset "ΔWoffdiag_element" begin

    mod = UEG()

    ens = CEnsemble(2.0, 5.680898543560106, 7)# θ: 0.125, λ: 0.09945178864947428

    t1 = ImgTime(0.3)
    t2 = ImgTime(0.7)
    τ3 = ImgTime(0.9)

    @test ΔWoffdiag_element(mod, ens, kinks_from_periodic_interval(sd, ImgTime(0.3), ImgTime(0.9)), kinks_from_periodic_interval(sd, ImgTime(0.3), ImgTime(0.7))) ≈ Woffdiag_element(mod, ens, d,c,b,a)
end


@testset "ΔWdiag_element(::Model, ::Ensemble, ::Configuration, i, j, k, l, τ1, τ2)" begin

    mod = UEG()

    i = PlaneWave((0,-4,0))
    j = PlaneWave((0,3,1))
    k = PlaneWave((0,0,1))
    l = PlaneWave((0,-1,0))

    τ1 = ImgTime(0.1)
    τ2 = ImgTime(0.3)
    τ3 = ImgTime(0.4)

    λ = 0.8
    β = 0.02
    N = length(conf.occupations)

    @test ΔWdiag_element(mod, CEnsemble(λ, β, N), conf, i, j, k, l, τ2, τ3) ≈ ΔW_diag(mod, i, j, k, l, occupations_at(conf,τ2)) * (τ3 - τ2) * λ

    @test ΔWdiag_element(mod, CEnsemble(λ, β, N), conf, i, j, k, l, τ1, τ3) ≈ ( ΔW_diag(mod, i, j, k, l, occupations_at(conf,τ1)) * (ImgTime(0.2) - τ1)
                                                                                + ΔW_diag(mod, i, j, k, l, occupations_at(conf,ImgTime(0.2))) * (τ3 - ImgTime(0.2))
                                                                                ) * λ

    @test ΔWdiag_element(mod, CEnsemble(λ, β, N), conf, i, j, k, l, τ3, τ1) ≈ ( ΔW_diag(mod, i, j, k, l, occupations_at(conf,τ3)) * (ImgTime(0.5) - τ3)
                                                                                + ΔW_diag(mod, i, j, k, l, occupations_at(conf,ImgTime(0.5))) * (ImgTime(0.6) - ImgTime(0.5))
                                                                                + ΔW_diag(mod, i, j, k, l, occupations_at(conf,ImgTime(0.6))) * (ImgTime(0.8) - ImgTime(0.6))
                                                                                + ΔW_diag(mod, i, j, k, l, occupations_at(conf,ImgTime(0.8))) * (ImgTime(1) + τ1 - ImgTime(0.8))
                                                                                ) * λ
end

@testset "τ_borders" begin
    @test τ_borders(conf, Set([a]), ImgTime(0.0)) == (ImgTime(0.8), ImgTime(0.2))
    @test τ_borders(conf, Set([a,e]), ImgTime(0.5)) == (ImgTime(0.2), ImgTime(0.6))
    @test τ_borders(conf, Set([b,d]), ImgTime(0.1)) == (ImgTime(0.8), ImgTime(0.2))
    @test τ_borders(conf, Set([e]), ImgTime(0.9)) == (ImgTime(0.0), ImgTime(1.0))

    conf1 = Configuration(sphere(PlaneWave((0,0,0),Up),dk=1))
    m = UEG()
    ens = CEnsemble(2,2,7)
    for _ in (1:5)
        dv, Δ = add_type_B(m, ens, conf1)
        if dv != 0
            apply_step!(conf1, Δ)
        end
        dv, Δ = add_type_C(m, ens, conf1)
        if dv != 0
            apply_step!(conf1, Δ)
        end
        dv, Δ = add_type_D(m, ens, conf1)
        if dv != 0
            apply_step!(conf1, Δ)
        end
        dv, Δ = add_type_E(m, ens, conf1)
        if dv != 0
            apply_step!(conf1, Δ)
        end
    end

    for _ in 1:10
        kink = rand(conf1.kinks)

        @test first(τ_borders(conf1, orbs(last(kink)), first(kink))) == τ_prev_affecting(conf1.kinks, orbs(last(kink)), first(kink))
        @test last(τ_borders(conf1, orbs(last(kink)), first(kink))) == τ_next_affecting(conf1.kinks, orbs(last(kink)), first(kink))
    end
end


@testset "add_orbs" begin
    c2 = Set{PlaneWave{3}}()
    @test add_orbs(c2,(a,)) == Set([a])
    @test add_orbs(c2,(a,b,c,)) == Set([a,b,c])
end

@testset "add_orbs!" begin
    c2 = Set{PlaneWave{3}}()
    add_orbs!(c2,(a,))
    @test c2 == Set([a])
    add_orbs!(c2,(b,c))
    @test c2 == Set([a,b,c])
end

@testset "drop_orbs" begin
    @test drop_orbs(Set([a,b,c]), (a,b)) == Set([c])
    @test drop_orbs(Set([a,b,c]), (c,)) == Set([a,b])
end

@testset "drop_orbs!" begin
    c2 = Set([a,b,c])
    drop_orbs!(c2, (a,))
    @test c2 == Set([b,c])
    drop_orbs!(c2, (b,c))
    @test c2 == Set{PlaneWave{3}}([])
end


@testset "add_kinks" begin
    @test add_kinks(Kinks( ImgTime(0.1) => T2(a,b), ImgTime(0.2) => T2(b,a), ImgTime(0.3) => T2(a,c) )
,(ImgTime(0.25) => T2(b,d),)) == Kinks( ImgTime(0.1) => T2(a,b), ImgTime(0.2) => T2(b,a), ImgTime(0.25) => T2(b,d), ImgTime(0.3) => T2(a,c) )
    @test add_kinks(Kinks( ImgTime(0.9) => T2(a,b), ImgTime(0.95) => T2(b,a) ), (ImgTime(0.1) => T2(a,b), ImgTime(0.2) => T2(b,a), ImgTime(0.3) => T2(a,c))) == Kinks( ImgTime(0.1) => T2(a,b), ImgTime(0.2) => T2(b,a), ImgTime(0.3) => T2(a,c), ImgTime(0.9) => T2(a,b), ImgTime(0.95) => T2(b,a) )
end

@testset "add_kinks!" begin
    c2 = Kinks( ImgTime(0.1) => T2(a,b), ImgTime(0.2) => T2(b,a), ImgTime(0.3) => T2(a,c) )
    add_kinks!(c2,(ImgTime(0.4) => T2(d,e),))
    @test c2 == Kinks( ImgTime(0.1) => T2(a,b), ImgTime(0.2) => T2(b,a), ImgTime(0.3) => T2(a,c), ImgTime(0.4) => T2(d,e))
    add_kinks!(c2,(ImgTime(0.9) => T2(a,b), ImgTime(0.95) => T2(b,a)))
    @test c2 == Kinks( ImgTime(0.1) => T2(a,b), ImgTime(0.2) => T2(b,a), ImgTime(0.3) => T2(a,c), ImgTime(0.4) => T2(d,e), ImgTime(0.9) => T2(a,b), ImgTime(0.95) => T2(b,a))
end

@testset "drop_kinks!" begin
    c2 = Kinks( ImgTime(0.1) => T2(a,b), ImgTime(0.2) => T2(b,a), ImgTime(0.3) => T2(a,c), ImgTime(0.4) => T2(d,e), ImgTime(0.9) => T2(a,b), ImgTime(0.95) => T2(b,a))
    drop_kinks!(c2, (ImgTime(0.1) => T2(a,b),))
    @test c2 == Kinks( ImgTime(0.2) => T2(b,a), ImgTime(0.3) => T2(a,c), ImgTime(0.4) => T2(d,e), ImgTime(0.9) => T2(a,b), ImgTime(0.95) => T2(b,a))
    drop_kinks!(c2, (ImgTime(0.4) => T2(d,e), ImgTime(0.9) => T2(a,b)))
    @test c2 == Kinks( ImgTime(0.2) => T2(b,a), ImgTime(0.3) => T2(a,c), ImgTime(0.95) => T2(b,a))
end

@testset "drop_kinks" begin
    c2 = Kinks( ImgTime(0.1) => T2(a,b), ImgTime(0.2) => T2(b,a), ImgTime(0.3) => T2(a,c), ImgTime(0.4) => T2(d,e), ImgTime(0.9) => T2(a,b), ImgTime(0.95) => T2(b,a))
    @test drop_kinks(c2, (ImgTime(0.1) => T2(a,b),)) == Kinks( ImgTime(0.2) => T2(b,a), ImgTime(0.3) => T2(a,c), ImgTime(0.4) => T2(d,e), ImgTime(0.9) => T2(a,b), ImgTime(0.95) => T2(b,a))
    @test drop_kinks(c2, (ImgTime(0.1) => T2(a,b), ImgTime(0.4) => T2(d,e), ImgTime(0.9) => T2(a,b))) == Kinks( ImgTime(0.2) => T2(b,a), ImgTime(0.3) => T2(a,c), ImgTime(0.95) => T2(b,a))
end

