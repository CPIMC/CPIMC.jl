
S = sphere_with_same_spin(OrbitalHEG((0,0,0)),dk=1)
a = OrbitalHEG((-2,0,0))
b = OrbitalHEG((3,0,0))
c = OrbitalHEG((0,0,0))
d = OrbitalHEG((1,0,0))
e = OrbitalHEG((5,9,9))

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


@testset "orbs_ordered(::T4)" begin
    @test orbs_ordered(T4(a,b,c,d))[1] == a
    @test orbs_ordered(T4(a,b,c,d))[2] == b
    @test orbs_ordered(T4(a,b,c,d))[3] == c
    @test orbs_ordered(T4(a,b,c,d))[4] == d
    @test orbs_ordered(T4(a,b,c,d)) == [a,b,c,d]
end

@testset "orbs_ordered(::SortedDict{ImgTime,<:Kink})" begin
    @test orbs_ordered(conf.kinks)[1] == a
    @test orbs_ordered(conf.kinks)[4] == d
    @test orbs_ordered(conf.kinks)[end] == a
    @test orbs_ordered(conf.kinks)[end-4] == c
end

#### tests for functions in `orbital.jl` #######################################
################################################################################

@testset "flip(o::OrbitalHEG)" begin
    @test flip(OrbitalHEG((1,0,0),Down)) == OrbitalHEG((1,0,0),Up)
    x = OrbitalHEG((0,),Down)
    @test flip(flip(x)) == x
end

@testset "dimension(o::OrbitalHEG{D}) where {D}" begin
    @test dimension(a) == 3
    @test dimension(OrbitalHEG((1,))) == 1
    @test dimension(conf.occupations) == 3
end

#### tests for formulas in `model.jl` ##########################################
################################################################################


@testset "fractional_spin_polarization(occ::Set{OrbitalHEG)" begin
    # equal number of Spin Up and Spin Down orbitals
    occ = Set([ OrbitalHEG((1,0,0),Down), OrbitalHEG((2,0,0),Up),
                OrbitalHEG((3,0,0),Down), OrbitalHEG((4,0,0),Up) ])
    @test fractional_spin_polarization(occ) == 0
    # 1 Spin Up and 3 Spin Down orbitals
    occ = Set([ OrbitalHEG((1,0,0),Down), OrbitalHEG((2,0,0),Down),
                OrbitalHEG((3,0,0),Down), OrbitalHEG((4,0,0),Up) ])
    @test fractional_spin_polarization(occ) == 0.5
    @test fractional_spin_polarization(occ) == fractional_spin_polarization(Set([flip(o) for o in occ]))
end

@testset "fractional_spin_polarization(c::Configuration{OrbitalHEG})" begin
    @test fractional_spin_polarization(conf) == 0
end

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
    @test β(Θ, r, N, ξ, d) ≈ (α * r)^2 * 16 / ( (1 + ξ)^(2/d) * (2π)^4 * λ(N, r, d)^2 * Θ )
end




@testset "compare parameter calculation with old versions for ξ=1 " begin

    " This is the expression used in the previous version. "
    function β_old(θ::Float64, N::Int, c::Union{Nothing,Configuration} = nothing )
        if isnothing(c) # no configuration provided
            Spin_Faktor = 1
        else
            Spin_Faktor = max(get_spin_up_down_count(c)...)/N
        end
        return ((2*pi)^2)/(((6*(pi^2)*N*Spin_Faktor)^(2/3))*θ)
    end

    "This is the expression as given in the thesis of T. Schoof."
    β_ref(θ::Float64, N::Int) = ((2*pi)^2)/(((6*(pi^2)*N)^(2/3))*θ)

    N = 18
    r = 0.56
    Θ = 0.123
    ξ = 1.0# use full polarization for comparison
    α = (4 / (9π))^(1/3)
    d = 3

    @test β_ref(Θ, N) ≈ β(Θ, r, N, ξ, d)

    " This is the expression used in the previous version. "
    function λ_old(N::Int, rs::Float64)
        return 4/((2*pi)^3) * (4*pi/3)^(1/3) * rs * N^(1/3)
    end

    @test λ_old(N, r) ≈ λ(N, r)

end
