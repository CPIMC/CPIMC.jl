include("../../src/UEG/model.jl")
e = Ensemble(2,β(0.125, 7),7)

S = sphere_with_same_spin(OrbitalHEG((0,0,0)),dk=1)
a = OrbitalHEG((-2,0,0))
b = OrbitalHEG((3,0,0))
c = OrbitalHEG((0,0,0))
d = OrbitalHEG((1,0,0))

sd = SortedDict{ImgTime, Kink{<:Orbital}}( ImgTime(0.5) => T4(a,b,c,d),
                                           ImgTime(0.6) => T4(c,d,a,b))

conf = Configuration(sphere(OrbitalHEG((0,0,0),Up),dk=1),sd)
conf_pol = Configuration(sphere_with_same_spin(OrbitalHEG((0,0,0),Up),dk=1),sd)



@testset "fractional_spin_polarization" begin
    @test fractional_spin_polarization(sphere_with_same_spin(OrbitalHEG((0,0,0),Up),dk=1)) == 1
    @test fractional_spin_polarization(sphere(OrbitalHEG((0,0,0),Up),dk=1)) == 0
    @test (β(0.125, 7, fractional_spin_polarization(sphere_with_same_spin(OrbitalHEG((0,0,0),Up),dk=1))) ==
        β(0.125, 14, fractional_spin_polarization(sphere(OrbitalHEG((0,0,0),Up),dk=1))))
end


@testset "wminus" begin
    r = 1:5

    for _ in 1:10000
        orb1 = OrbitalHEG((rand(r),rand(r),rand(r)),rand([Up,Down]))
        orb2 = orb1
        while orb1 == orb2
            orb2 = OrbitalHEG((rand(r),rand(r),rand(r)),rand([Up,Down]))
        end
        orb3 = orb1
        while (orb3 == orb2) || (orb3 == orb1)
            orb3 = OrbitalHEG((rand(r),rand(r),rand(r)),rand([orb1.spin,orb2.spin]))
        end

        if orb1.spin == orb2.spin
            orb4 = OrbitalHEG((orb1.vec + orb2.vec - orb3.vec), orb3.spin)
        else
            orb4 = OrbitalHEG((orb1.vec + orb2.vec - orb3.vec), flip(orb3.spin))
        end

        @test !isinf(abs(wminus(orb1,orb2,orb3,orb4)))
        @test !isnan(abs(wminus(orb1,orb2,orb3,orb4)))
    end
end


@testset "W_diag" begin
    τ1 = ImgTime(0.8)
    τ2 = ImgTime(0.2)
    Δ = Step(Set{OrbitalHEG{3}}([c,d]), Configuration(Set{OrbitalHEG{3}}([a,b]), τ1 => T4(b,a,c,d), τ2 => T4(d,c,b,a)))
    @test round(W_diag(e,apply_step(conf_pol, Δ)) - W_diag(e,conf_pol),digits = 11) == round(Δdiagonal_interaction(conf_pol, e, a::Orbital, b::Orbital, c::Orbital, d::Orbital, τ1, τ2),digits= 11)
    @test round(W_diag(e,apply_step(conf, Δ)) - W_diag(e,conf),digits = 11) == round(Δdiagonal_interaction(conf, e, a::Orbital, b::Orbital, c::Orbital, d::Orbital, τ1, τ2), digits= 11)
end
