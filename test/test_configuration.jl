
S = get_sphere_with_same_spin(OrbitalHEG((0,0,0),1),dk=1)
a = OrbitalHEG((-2,0,0),1)
b = OrbitalHEG((3,0,0),1)
c = OrbitalHEG((0,0,0),1)
d = OrbitalHEG((1,0,0),1)
e = OrbitalHEG((5,9,9),1)

sd = SortedDict{ImgTime, Kink{<:Orbital}}( ImgTime(0.2) => T4(a,b,c,d),
                                           ImgTime(0.5) => T4(c,d,a,b),
                                           ImgTime(0.6) => T4(b,a,d,c),
                                           ImgTime(0.8) => T4(d,c,b,a) )

conf = Configuration(get_sphere(OrbitalHEG((0,0,0),1),dk=1),sd)

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


@testset "orbs_ordered(::T4)" begin
    @test orbs_ordered(T4(a,b,c,d))[1] == a
    @test orbs_ordered(T4(a,b,c,d))[2] == b
    @test orbs_ordered(T4(a,b,c,d))[3] == c
    @test orbs_ordered(T4(a,b,c,d))[4] == d
    @test orbs_ordered(T4(a,b,c,d)) == [a,b,c,d]
end

@testset "orbs_ordered(::SortedDict{ImgTime,<:Kink})" begin
    orbs_ordered(conf.kinks)[1] == a
    orbs_ordered(conf.kinks)[4] == d
    orbs_ordered(conf.kinks)[end] == a
    orbs_ordered(conf.kinks)[end-4] == d
end

println(typeof(orbs_ordered(conf.kinks)))
