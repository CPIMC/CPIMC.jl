
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
    @test β(Θ, N, ξ, d) ≈ (α * r)^2 * 16 / ( (1 + ξ)^(2/d) * (2π)^4 * λ(N, r, d)^2 * Θ )
end

@testset "compare parameter calculation with old versions " begin

    " returns a tuple of number of particles with spin up and number of particles with spin down "
    function get_spin_up_down_count(c::Configuration)
        up = 0
        down = 0
        for orb in c.occupations
            if orb.spin == Up
                up += 1
            elseif orb.spin == Down
                down += 1
            else
                @assert(false)#Spin is not 1 or -1
            end
        end
        return (up,down)
    end

    " This is the expression used in the previous version. "
    function β_old(θ::Float64, N::Int, c::Union{Nothing,Configuration} = nothing )
        if isnothing(c) # no configuration provided
            Spin_Faktor = 1
        else
            Spin_Faktor = max(get_spin_up_down_count(c)...)/N
        end
        #println("Spin_Faktor : ", Spin_Faktor)
        return ((2*pi)^2)/(((6*(pi^2)*N*Spin_Faktor)^(2/3))*θ)
    end

    "This is the expression as given in the thesis of T. Schoof."
    β_ref(θ::Float64, N::Int) = ((2*pi)^2)/(((6*(pi^2)*N)^(2/3))*θ)

    N = 4
    r = 0.56
    Θ = 0.123
    ξ = 1.0# use full polarization for comparison
    α = (4 / (9π))^(1/3)
    d = 3

    @test β_ref(Θ, N) ≈ β(Θ, N, ξ, d)

    ## test for unpolarized system
    c = Configuration(Set([OrbitalHEG((1,),Down), OrbitalHEG((2,),Up), OrbitalHEG((3,),Down), OrbitalHEG((4,),Up)]))


    #println("spin_count : ", get_spin_up_down_count(c))
    #println("xi = ", fractional_spin_polarization(c))

    @test β_old(Θ, N, c) ≈ β(Θ, N, fractional_spin_polarization(c), d)

    " This is the expression used in the previous version. "
    function λ_old(N::Int, rs::Float64)
        return 4/((2*pi)^3) * (4*pi/3)^(1/3) * rs * N^(1/3)
    end

    @test λ_old(N, r) ≈ λ(N, r)

    """
        β(θ::Float64, N::Int, ξ::Float64 )

    calculate β in internal units
    """
    function β_previus_version(θ::Float64, N::Int, ξ::Float64 = 1.0)
        return (2*pi)^2/(((6*(pi^2) * N/2 * (1+abs(ξ)))^(2/3))*θ)
    end


    for i in 1:100
        θ = rand()*100
        rs = rand()*100
        N = rand(1:100)
        ξ = 1.0
        for i in 1:10
            #@test β(θ, N, ξ, rs) ≈ β(θ, N, ξ, rs*rand()*100)
            #@test β(θ, N, ξ, rs,2) ≈ β(θ, N, ξ, rs*rand()*100,2)
        end
        @test (β_old(θ, N) ≈ β(θ, N, ξ))
        @test (β_ref(θ, N) ≈ β(θ, N, ξ))
        ξ = rand()
        @test (β_previus_version(θ, N, ξ) ≈ β(θ, N, ξ))
    end

end


e = Ensemble(2,5.68089,7)

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
