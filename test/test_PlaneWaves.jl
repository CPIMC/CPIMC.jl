using CPIMC, DataStructures
import CPIMC: orbs, adjacent_kinks_affecting_orbs, kinks_affecting_orbs, τ_borders, isunaffected, time_ordered_orbs, occupations_at, longest_type_1_chain_length, right_type_1_count
using CPIMC.PlaneWaves

a = PlaneWave((-2,0,0))
b = PlaneWave((3,0,0))
c = PlaneWave((0,0,0))
d = PlaneWave((1,0,0))
e = PlaneWave((5,9,9))

sd = Kinks( ImgTime(0.2) => T4(a,b,c,d),
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
    orb4 =  find_fourth_orb_for_kink(orb3, orb1, orb2)
    @test find_fourth_orb_for_kink(orb3, orb1, orb2) == find_fourth_orb_for_kink(orb3, orb3, orb4)
end
