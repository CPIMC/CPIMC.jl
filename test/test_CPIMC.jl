using CPIMC, CPIMC.PlaneWaves, CPIMC.UniformElectronGas, CPIMC.DefaultUpdates, DataStructures

import CPIMC: orbs, adjacent_kinks_affecting_orbs, kinks_affecting_orbs, prev, next, next_affecting, prev_affecting, τ_prev_affecting, τ_next_affecting, τ_borders, isunaffected, isunaffected_in_interval, time_ordered_orbs, occupations_at, longest_type_1_chain_length, right_type_1_count, kinks_from_periodic_interval, times_from_periodic_interval, Δ, Woffdiag_element, ΔWoffdiag_element, ΔWdiag_element, ΔW_diag, add_orbs, add_orbs!, drop_orbs, drop_orbs!, drop_kinks, drop_kinks!, add_kinks, add_kinks!

S = sphere_with_same_spin(PlaneWave((0,0,0)),dk=1)
a = PlaneWave((-2,0,0))
b = PlaneWave((3,0,0))
c = PlaneWave((0,0,0))
d = PlaneWave((1,0,0))
e = PlaneWave((5,9,9))
f = PlaneWave((1,1,1))

g = PlaneWave(a.vec + e.vec - f.vec, Up)
h = PlaneWave(c.vec + d.vec - g.vec, Up)

sd = Kinks( ImgTime(0.2) => T4(a,b,c,d),
       ImgTime(0.5) => T4(c,d,a,b),
       ImgTime(0.6) => T4(b,a,d,c),
       ImgTime(0.8) => T4(d,c,b,a) )

conf = Configuration(S,sd)

@testset "apply_step" begin
    @test apply_step(conf, Step()) == conf
    @test apply_step(conf, Step((c,d),nothing,(e,f),nothing)).occupations == union(setdiff(S,(c,d)),(e,f))
end

@testset "apply_step!" begin
    c2 = Configuration(S,sd)
    apply_step!(c2, Step())
    @test c2 == conf
    apply_step!(c2, Step((c,d),nothing,(e,f),nothing))
    @test c2.occupations == union(setdiff(S,(c,d)),(e,f))
end




@testset "apply_step" for _ in (1:200)
    c1 = Configuration(sphere(PlaneWave((0,0,0),Up),dk=1))
    c2 = Configuration(sphere(PlaneWave((0,0,0),Up),dk=1))
    m = UEG()
    ens = CEnsemble(2,2,7)
    step_list = Array{Step,1}()
    for _ in (1:5)
        dv, Δ = add_type_B(m, ens, c1)

        if dv != 0
            push!(step_list, Δ)
            apply_step!(c1, Δ)
        end
        dv, Δ = add_type_C(m, ens, c1)

        if dv != 0
            push!(step_list, Δ)
            apply_step!(c1, Δ)
        end
        dv, Δ = add_type_D(m, ens, c1)

        if dv != 0
            push!(step_list, Δ)
            apply_step!(c1, Δ)
        end
        dv, Δ = add_type_E(m, ens, c1)

        if dv != 0
            push!(step_list, Δ)
            apply_step!(c1, Δ)
        end
    end
    for i in (1:3)
        dv, Δ = remove_type_B(m, ens, c1)

        if dv != 0
            push!(step_list, Δ)
            apply_step!(c1, Δ)
        end
        dv, Δ = remove_type_C(m, ens, c1)

        if dv != 0
                push!(step_list, Δ)
            apply_step!(c1, Δ)
        end
        dv, Δ = remove_type_D(m, ens, c1)

        if dv != 0
            push!(step_list, Δ)
            apply_step!(c1, Δ)
        end
        dv, Δ = remove_type_E(m, ens, c1)

        if dv != 0
            push!(step_list, Δ)
            apply_step!(c1, Δ)
        end
    end
    @test apply_step(c2,step_list).kinks == c1.kinks
    @test apply_step(c2,step_list).occupations == c1.occupations
end
