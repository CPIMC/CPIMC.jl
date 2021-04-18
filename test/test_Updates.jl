using CPIMC, CPIMC.PlaneWaves, DataStructures, CPIMC.DefaultUpdates
import LinearAlgebra: dot
import CPIMC: orbs, adjacent_kinks_affecting_orbs, kinks_affecting_orbs, τ_borders, isunaffected, time_ordered_orbs, occupations_at, longest_type_1_chain_length, right_type_1_count

@testset "removable_pairs" for _ in (1:50)
    conf1 = Configuration(sphere(PlaneWave((0,0,0),Up),dk=1))
    m = UEG()
    ens = CEnsemble(2,2,7)
    for _ in (1:5)
        dv, Δ = add_type_B(m, ens, conf1)
        if dv != 0
            apply_step!(conf1, Δ)
        end
        dv, Δ = add_type_C(m, ens, conf1)
        if dv != 0
            apply_step!(conf1, Δ)
        end
        dv, Δ = add_type_D(m, ens, conf1)
        if dv != 0
            apply_step!(conf1, Δ)
        end
        dv, Δ = add_type_E(m, ens, conf1)
        if dv != 0
            apply_step!(conf1, Δ)
        end
    end
    for pair in CPIMC.DefaultUpdates.right_type_B_removable_pairs(conf1.kinks)
        (_,kink1) = conf1.kinks[first(pair)]
        @test dot(kink1.i.vec-kink1.k.vec, kink1.i.vec-kink1.k.vec) <= ex_radius^2
        @test kink1.i.spin == kink1.k.spin
    end
    for pair in CPIMC.DefaultUpdates.right_type_C_removable_pairs(conf1.kinks)
        (_,kink1) = conf1.kinks[first(pair)]
        @test dot(kink1.i.vec-kink1.k.vec, kink1.i.vec-kink1.k.vec) <= ex_radius^2
        @test kink1.i.spin == kink1.k.spin
    end
    for pair in CPIMC.DefaultUpdates.left_type_C_removable_pairs(conf1.kinks)
        (_,kink1) = conf1.kinks[first(pair)]
        @test dot(kink1.i.vec-kink1.k.vec, kink1.i.vec-kink1.k.vec) <= ex_radius^2
        @test kink1.i.spin == kink1.k.spin
    end
    for pair in CPIMC.DefaultUpdates.right_type_D_removable_pairs(conf1.kinks)
        (_,kink1) = conf1.kinks[first(pair)]
        @test dot(kink1.i.vec-kink1.k.vec, kink1.i.vec-kink1.k.vec) <= ex_radius^2
        @test kink1.i.spin == kink1.k.spin
    end
    for pair in CPIMC.DefaultUpdates.left_type_D_removable_pairs(conf1.kinks)
        (_,kink1) = conf1.kinks[first(pair)]
        @test dot(kink1.i.vec-kink1.k.vec, kink1.i.vec-kink1.k.vec) <= ex_radius^2
        @test kink1.i.spin == kink1.k.spin
    end
    for pair in CPIMC.DefaultUpdates.right_type_E_removable_pairs(conf1.kinks)
        kink1 = last(first(pair))
        @test dot(kink1.i.vec-kink1.k.vec, kink1.i.vec-kink1.k.vec) <= ex_radius^2
    end
    for pair in CPIMC.DefaultUpdates.left_type_E_removable_pairs(conf1.kinks)
        kink1 = last(first(pair))
        @test dot(kink1.i.vec-kink1.k.vec, kink1.i.vec-kink1.k.vec) <= ex_radius^2
    end

end


@testset "possible_new_orb" for _ in (1:50)
    conf = Configuration(sphere(PlaneWave((0,0,0),Up),dk=1))
    m = UEG()
    ens = CEnsemble(2,2,7)
    for _ in (1:5)
        dv, Δ = add_type_B(m, ens, conf)
        if dv != 0
            apply_step!(conf, Δ)
        end
        dv, Δ = add_type_C(m, ens, conf)
        if dv != 0
            apply_step!(conf, Δ)
        end
        dv, Δ = add_type_D(m, ens, conf)
        if dv != 0
            apply_step!(conf, Δ)
        end
        dv, Δ = add_type_E(m, ens, conf)
        if dv != 0
            apply_step!(conf, Δ)
        end
    end
    old_kink = rand(conf.kinks)
    occs = CPIMC.occupations_at(conf, first(old_kink))


    opportunities_new_orb1 = CPIMC.DefaultUpdates.possible_new_orb1_C(occs, last(old_kink).k, last(old_kink).l,last(old_kink).i, last(old_kink).j)
    for new_orb1 in opportunities_new_orb1
        new_orb2 = find_fourth_orb_for_kink(new_orb1, last(old_kink).i, last(old_kink).j)
        @test dot(new_orb1.vec-last(old_kink).k.vec, new_orb1.vec-last(old_kink).k.vec) <= ex_radius^2
        @test in(new_orb1.spin, [last(old_kink).i.spin, last(old_kink).j.spin])
        @test (new_orb1.spin == new_orb2.spin) == (last(old_kink).i.spin == last(old_kink).j.spin)
        @test(new_orb1.vec + new_orb2.vec == last(old_kink).i.vec + last(old_kink).j.vec)
        @test length(Set([new_orb2, new_orb1, last(old_kink).i, last(old_kink).j])) == 4
    end

    opportunities_new_orb1 = CPIMC.DefaultUpdates.possible_new_orb1_D(occs, last(old_kink).k, last(old_kink).l,last(old_kink).i, last(old_kink).j)
    for new_orb1 in opportunities_new_orb1
        new_orb2 = find_fourth_orb_for_kink(new_orb1, last(old_kink).i, last(old_kink).j)
        @test dot(new_orb1.vec-last(old_kink).k.vec, new_orb1.vec-last(old_kink).k.vec) <= ex_radius^2
        @test in(new_orb1.spin, [last(old_kink).i.spin, last(old_kink).j.spin])
        @test (new_orb1.spin == new_orb2.spin) == (last(old_kink).i.spin == last(old_kink).j.spin)
        @test(new_orb1.vec + new_orb2.vec == last(old_kink).i.vec + last(old_kink).j.vec)
        @test length(Set([new_orb2, new_orb1, last(old_kink).i, last(old_kink).j])) == 4
    end

    for new_kink_old_creator in [last(old_kink).i, last(old_kink).j]
        changed_kink_old_creator = setdiff([last(old_kink).i, last(old_kink).j], [new_kink_old_creator])[1]
        for new_kink_old_annihilator in [last(old_kink).k, last(old_kink).l]
            changed_kink_old_annihilator = setdiff([last(old_kink).k, last(old_kink).l], [new_kink_old_annihilator])[1]

            opportunities_new_kink_new_annihilator = CPIMC.DefaultUpdates.possible_new_kink_new_occ_orb(occs, new_kink_old_creator, changed_kink_old_creator, new_kink_old_annihilator, changed_kink_old_annihilator)
            for new_kink_new_annihilator in opportunities_new_kink_new_annihilator
                new_kink_new_creator = find_fourth_orb_for_kink(new_kink_old_creator, new_kink_new_annihilator, new_kink_old_annihilator)
                @test dot(new_kink_new_annihilator.vec-new_kink_old_creator.vec, new_kink_new_annihilator.vec-new_kink_old_creator.vec) <= ex_radius^2
                @test in(new_kink_new_annihilator.spin, [new_kink_old_creator.spin, new_kink_new_creator.spin])
                @test (new_kink_new_annihilator.spin == new_kink_old_annihilator.spin) == (new_kink_old_creator.spin == new_kink_new_creator.spin)
                @test(new_kink_new_annihilator.vec + new_kink_old_annihilator.vec == new_kink_old_creator.vec + new_kink_new_creator.vec)
                @test length(Set([new_kink_new_creator, new_kink_new_annihilator, new_kink_old_creator, new_kink_old_annihilator])) == 4
            end
        end
    end
end
