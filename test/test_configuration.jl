using CPIMC, CPIMC.PlaneWaves, DataStructures
import CPIMC: ImgTime, orbs, T2, T4, adjacent_kinks_affecting_orbs, kinks_affecting_orbs, τ_borders, isunaffected, time_ordered_orbs, occupations, longest_type_1_chain_length, right_type_1_count


S = sphere_with_same_spin(PlaneWave((0,0,0)),dk=1)
a = PlaneWave((-2,0,0))
b = PlaneWave((3,0,0))
c = PlaneWave((0,0,0))
d = PlaneWave((1,0,0))
e = PlaneWave((5,9,9))
f = PlaneWave((1,1,1))

g = PlaneWave(a.vec + e.vec - f.vec, Up)
h = PlaneWave(c.vec + d.vec - g.vec, Up)



sd = SortedDict{ImgTime, Kink{<:Orbital}}( ImgTime(0.2) => T4(a,b,c,d),
                                           ImgTime(0.5) => T4(c,d,a,b),
                                           ImgTime(0.6) => T4(b,a,d,c),
                                           ImgTime(0.8) => T4(d,c,b,a) )

conf = Configuration(sphere(PlaneWave((0,0,0),Up),dk=1),sd)

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

@testset "Type_1_investigation" begin
    g = PlaneWave(a.vec + e.vec - f.vec, Up)
    h = PlaneWave(c.vec + d.vec - g.vec, Up)

    Type_1_chain = SortedDict{ImgTime, Kink{<:Orbital}}( ImgTime(0.2) => T4(a,b,c,d),
                                               ImgTime(0.5) => T4(f,g,e,a),
                                               ImgTime(0.6) => T4(c,d,h,g),
                                               ImgTime(0.8) => T4(e,h,b,f) )

    @test (a.vec + b.vec - c.vec - d.vec) == PlaneWave((0,0,0)).vec
    @test (f.vec + g.vec - e.vec - a.vec) == PlaneWave((0,0,0)).vec
    @test (c.vec + d.vec - h.vec - g.vec) == PlaneWave((0,0,0)).vec
    @test (e.vec + h.vec - b.vec - f.vec) == PlaneWave((0,0,0)).vec

    occs = setdiff!(union!(sphere(PlaneWave((0,0,0),Up),dk=1),
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
