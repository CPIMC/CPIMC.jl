using CPIMC
using CPIMC.UniformElectronGas
import CPIMC.UniformElectronGas: w

import CPIMC: orbs

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






@testset "ΔW_diag for polarized occupation" begin

    mod = UEG()

    a = PlaneWave((-2,0,0))
    b = PlaneWave((3,0,0))
    c = PlaneWave((0,0,0))
    d = PlaneWave((1,0,0))

    ## occupation to be changed
    occ = Set([a,b,c,d])

    ### Test 2-particle excitation
    # choose creator orbitals
    i = PlaneWave((1,1,1))
    j = PlaneWave((0,-1,-1))
    # choose annihilator orbitals
    k = c
    l = d
    @assert iszero( i.vec + j.vec - k.vec - l.vec ) " momentum not conserved for this excitation "
    @assert (i.spin == k.spin) & (j.spin == l.spin) " spin is not conserved for this excitation "

    # this functions sums all combinations of orbitals o ∈ occ with annihilated orbitals k (o ≠ k) and l (o ≠ l)
    # and substracts the sum of all combinations of orbitals o ∈ ( occ ̸ {k,l} ) ∪ {i,j} of the new occupation created with orbitals i (o ≠ i) and j (o ≠ j)
    @test ΔW_diag(mod, i, j, k, l, occ) ≈ ( w(mod, a,k,k,a) + w(mod, b,k,k,b) + w(mod, a,l,l,a) + w(mod, b,l,l,b) + w(mod, k,l,l,k)
                                      - w(mod, a,i,i,a) - w(mod, b,i,i,b) - w(mod, a,j,j,a) - w(mod, b,j,j,b) - w(mod, i,j,j,i) )

    ### Test 1-particle excitation
    i = PlaneWave((-3,0,0))
    j = b
    # this functions sums all combinations of orbitals o ∈ occ with the annihilated orbital j (o ≠ j)
    # and substracts the sum of all combinations of orbitals o ∈ ( occ ̸ {j} ) ∪ {i} of the new occupation created with orbital i (o ≠ i)
    @test ΔW_diag(mod, i, j, occ) ≈ ( w(mod, a,j,j,a) + w(mod, c,j,j,c) + w(mod, d,j,j,d)
                                - w(mod, a,i,i,a) - w(mod, c,i,i,c) - w(mod, d,i,i,d) )
end


@testset "ΔW_diag for unpolarized occupation" begin

    mod = UEG()

    a = PlaneWave((-2,0,0),Up)
    b = PlaneWave((3,0,0),Down)
    c = PlaneWave((0,0,0),Up)
    d = PlaneWave((1,0,0),Down)

    ## occupation to be changed
    occ = Set([a,b,c,d])

    ### Test 2-particle excitation
    # choose creator orbitals
    i = PlaneWave((1,1,1),Up)
    j = PlaneWave((0,-1,-1),Down)
    # choose annihilator orbitals
    k = c
    l = d
    @assert iszero( i.vec + j.vec - k.vec - l.vec ) " momentum not conserved for this excitation "
    @assert (i.spin == k.spin) & (j.spin == l.spin) " spin is not conserved for this excitation "

    # this functions sums all combinations of orbitals o ∈ occ with annihilated orbitals k (o ≠ k) and l (o ≠ l)
    # and substracts the sum of all combinations of orbitals o ∈ ( occ ̸ {k,l} ) ∪ {i,j} of the new occupation created with orbitals i (o ≠ i) and j (o ≠ j)
    @test ΔW_diag(mod, i, j, k, l, occ) ≈ ( w(mod, a,k,k,a) + w(mod, b,k,k,b) + w(mod, a,l,l,a) + w(mod, b,l,l,b) + w(mod, k,l,l,k)
                                      - w(mod, a,i,i,a) - w(mod, b,i,i,b) - w(mod, a,j,j,a) - w(mod, b,j,j,b) - w(mod, i,j,j,i) )

    ## test for orbitals where a momentum vector of one creator orbital is equal to the momentum vector of one orbital in the occupation
    # choose creator orbitals
    i = PlaneWave((-2,0,0),Down)
    j = PlaneWave((5,0,0),Up)
    # choose annihilator orbitals
    k = b
    l = c
    @assert iszero( i.vec + j.vec - k.vec - l.vec ) " momentum not conserved for this excitation "
    @assert (i.spin == k.spin) & (j.spin == l.spin) " spin is not conserved for this excitation "

    # the contribution w(a,i,i,a) does not occur here because a.vec == i.vec
    @test ΔW_diag(mod, i, j, k, l, occ) ≈ ( w(mod, a,k,k,a) + w(mod, d,k,k,d) + w(mod, a,l,l,a) + w(mod, d,l,l,d) + w(mod, k,l,l,k)
                                      # - w(a,i,i,a) does not occur here because a.vec == i.vec
                                       - w(mod, d,i,i,d) - w(mod, a,j,j,a) - w(mod, d,j,j,d) - w(mod, i,j,j,i) )

    ### Test 1-particle excitation
    i = PlaneWave((-3,0,0),Down)
    j = b
    # this functions sums all combinations of orbitals o ∈ occ with the annihilated orbital j (o ≠ j)
    # and substracts the sum of all combinations of orbitals o ∈ ( occ ̸ {j} ) ∪ {i} of the new occupation created with orbital i (o ≠ i)
    @test ΔW_diag(mod, i, j, occ) ≈ ( w(mod, a,j,j,a) + w(mod, c,j,j,c) + w(mod, d,j,j,d)
                                - w(mod, a,i,i,a) - w(mod, c,i,i,c) - w(mod, d,i,i,d) )
end
