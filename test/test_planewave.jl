using CPIMC, DataStructures
import CPIMC: ImgTime, orbs, T2, T4, adjacent_kinks_affecting_orbs, kinks_affecting_orbs, τ_borders, isunaffected, time_ordered_orbs, occupations, longest_type_1_chain_length, right_type_1_count

a = PlaneWave((-2,0,0))
b = PlaneWave((3,0,0))
c = PlaneWave((0,0,0))
d = PlaneWave((1,0,0))
e = PlaneWave((5,9,9))

sd = SortedDict{ImgTime, Kink{<:Orbital}}( ImgTime(0.2) => T4(a,b,c,d),
                                           ImgTime(0.5) => T4(c,d,a,b),
                                           ImgTime(0.6) => T4(b,a,d,c),
                                           ImgTime(0.8) => T4(d,c,b,a) )

conf = Configuration(sphere(PlaneWave((0,0,0),Up),dk=1),sd)


@testset "flip(o::PlaneWave)" begin
    @test flip(PlaneWave((1,0,0),Down)) == PlaneWave((1,0,0),Up)
    x = PlaneWave((0,),Down)
    @test flip(flip(x)) == x
end

@testset "dimension(o::PlaneWave{D}) where {D}" begin
    @test dimension(a) == 3
    @test dimension(PlaneWave((1,))) == 1
    @test dimension(conf.occupations) == 3
end

@testset "fractional_spin_polarization(occ::Set{PlaneWave)" begin
    # equal number of Spin Up and Spin Down orbitals
    occ = Set([ PlaneWave((1,0,0),Down), PlaneWave((2,0,0),Up),
                PlaneWave((3,0,0),Down), PlaneWave((4,0,0),Up) ])
    @test fractional_spin_polarization(occ) == 0
    # 1 Spin Up and 3 Spin Down orbitals
    occ = Set([ PlaneWave((1,0,0),Down), PlaneWave((2,0,0),Down),
                PlaneWave((3,0,0),Down), PlaneWave((4,0,0),Up) ])
    @test fractional_spin_polarization(occ) == 0.5
    @test fractional_spin_polarization(occ) == fractional_spin_polarization(Set([flip(o) for o in occ]))
end

@testset "fractional_spin_polarization(c::Configuration{PlaneWave})" begin
    @test fractional_spin_polarization(conf) == 0
end


