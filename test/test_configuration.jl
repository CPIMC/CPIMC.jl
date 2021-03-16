include("../src/Configuration.jl")
include("../src/UEG/model.jl")

using Test


S = get_sphere_with_same_spin(OrbitalHEG((0,0,0),1),dk=1)
a = OrbitalHEG((-2,0,0),1)
b = OrbitalHEG((3,0,0),1)
c = OrbitalHEG((0,0,0),1)
d = OrbitalHEG((1,0,0),1)

conf = Configuration(get_sphere(OrbitalHEG((0,0,0),1),dk=1),
        ImgTime(0.2) => T4(a,b,c,d),
        ImgTime(0.5) => T4(c,d,a,b),
        ImgTime(0.6) => T4(b,a,d,c),
        ImgTime(0.8) => T4(d,c,b,a))

@testset "adjacent_kinks_affecting_orbs" begin
    @test adjacent_kinks_affecting_orbs(conf, Set([a]), ImgTime(0.1)) == (ImgTime(0.8) => T4(d,c,b,a), ImgTime(0.2) => T4(a,b,c,d))
end
