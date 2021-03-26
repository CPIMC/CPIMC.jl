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
        println("Spin_Faktor : ", Spin_Faktor)
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

    @test β_ref(Θ, N) ≈ β(Θ, r, N, ξ, d)

    ## test for unpolarized system
    c = Configuration(Set([OrbitalHEG((1,),Down), OrbitalHEG((2,),Up), OrbitalHEG((3,),Down), OrbitalHEG((4,),Up)]))


    println("spin_count : ", get_spin_up_down_count(c))
    println("xi = ", fractional_spin_polarization(c))

    @test β_old(Θ, N, c) ≈ β(Θ, r, N, fractional_spin_polarization(c), d)

    " This is the expression used in the previous version. "
    function λ_old(N::Int, rs::Float64)
        return 4/((2*pi)^3) * (4*pi/3)^(1/3) * rs * N^(1/3)
    end

    @test λ_old(N, r) ≈ λ(N, r)

end

