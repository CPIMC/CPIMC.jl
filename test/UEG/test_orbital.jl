
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

