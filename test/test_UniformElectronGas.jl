using CPIMC
using CPIMC.UniformElectronGas
using CPIMC.Estimators

import CPIMC: ImgTime, orbs, T2, T4, Step, wdiag, wminus, apply_step, Δdiagonal_interaction

@testset "λ(N, rs, d) and rs(N::Int, λ::Float64, d) are consistent" begin
    N = 18
    r = 0.56
    @test rs(N, λ(N, r, 3), 3) ≈ r
    @test rs(N, λ(N, r, 2), 2) ≈ r
end


@testset " consistence of β(Θ::Float64, rs::Float64, ξ, d) with λ(N, rs, d) and rs(N::Int, λ::Float64, d) " begin
    N = 18
    r = 0.56
    Θ = 0.123
    ξ = 0.254
    α = (4 / (9π))^(1/3)
    d = 3
    @test β(Θ, N, ξ, d) ≈ (α * r)^2 * 16 / ( (1 + ξ)^(2/d) * (2π)^4 * λ(N, r, d)^2 * Θ )
end



e = CEnsemble(2,5.68089,7)

S = sphere_with_same_spin(PlaneWave((0,0,0)),dk=1)
a = PlaneWave((-2,0,0))
b = PlaneWave((3,0,0))
c = PlaneWave((0,0,0))
d = PlaneWave((1,0,0))

sd = SortedDict{ImgTime, Kink{<:Orbital}}( ImgTime(0.5) => T4(a,b,c,d),
                                           ImgTime(0.6) => T4(c,d,a,b))

conf = Configuration(sphere(PlaneWave((0,0,0),Up),dk=1),sd)
conf_pol = Configuration(sphere_with_same_spin(PlaneWave((0,0,0),Up),dk=1),sd)

m = UEG()

@testset "wminus" begin
    r = 1:5

    for _ in 1:10000
        orb1 = PlaneWave((rand(r),rand(r),rand(r)),rand([Up,Down]))
        orb2 = orb1
        while orb1 == orb2
            orb2 = PlaneWave((rand(r),rand(r),rand(r)),rand([Up,Down]))
        end
        orb3 = orb1
        while (orb3 == orb2) || (orb3 == orb1)
            orb3 = PlaneWave((rand(r),rand(r),rand(r)),rand([orb1.spin,orb2.spin]))
        end

        if orb1.spin == orb2.spin
            orb4 = PlaneWave((orb1.vec + orb2.vec - orb3.vec), orb3.spin)
        else
            orb4 = PlaneWave((orb1.vec + orb2.vec - orb3.vec), flip(orb3.spin))
        end

        @test !isinf(abs(wminus(m, orb1, orb2, orb3, orb4)))
        @test !isnan(abs(wminus(m, orb1, orb2, orb3, orb4)))
    end
end

@testset "wdiag" begin
    τ1 = ImgTime(0.8)
    τ2 = ImgTime(0.2)
    Δ = Step(Set{PlaneWave{3}}([c,d]), Configuration(Set{PlaneWave{3}}([a,b]), τ1 => T4(b,a,c,d), τ2 => T4(d,c,b,a)))
    @test round(W_diag(m, e, apply_step(conf_pol, Δ)) - W_diag(m, e,conf_pol),digits = 11) == round(Δdiagonal_interaction(m, e, conf_pol, a::Orbital, b::Orbital, c::Orbital, d::Orbital, τ1, τ2),digits= 11)
    @test round(W_diag(m, e, apply_step(conf, Δ)) - W_diag(m, e, conf),digits = 11) == round(Δdiagonal_interaction(m, e, conf, a::Orbital, b::Orbital, c::Orbital, d::Orbital, τ1, τ2), digits= 11)
end
