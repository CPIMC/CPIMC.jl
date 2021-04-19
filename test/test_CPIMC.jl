using CPIMC, CPIMC.PlaneWaves, CPIMC.UniformElectronGas, CPIMC.DefaultUpdates, DataStructures
import CPIMC: orbs, adjacent_kinks_affecting_orbs, kinks_affecting_orbs, prev, next, prev_index, next_index, index_prev_affecting, index_next_affecting, next_affecting, prev_affecting, τ_prev_affecting, τ_next_affecting, τ_borders, isunaffected, isunaffected_in_interval, time_ordered_orbs, occupations_at, longest_type_1_chain_length, right_type_1_count, kinks_from_periodic_interval, times_from_periodic_interval, Δ, Woffdiag_element, ΔWoffdiag_element, ΔWdiag_element, ΔW_diag


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
    for i in (1:2)
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


#=
#@testset "apply_step" for _ in (1:1)
function test()
    c1 = Configuration(sphere(PlaneWave((0,0,0),Up),dk=1))
    c2 = Configuration(sphere(PlaneWave((0,0,0),Up),dk=1))
    m = UEG()
    ens = CEnsemble(2,2,7)
    step_list = Array{Step,1}()
    for _ in (1:1)
        dv, Δ = add_type_B(m, ens, c1)
        push!(step_list, Δ)
        apply_step!(c1, Δ)
        dv, Δ = add_type_C(m, ens, c1)
        push!(step_list, Δ)
        apply_step!(c1, Δ)
        dv, Δ = add_type_D(m, ens, c1)
        push!(step_list, Δ)
        apply_step!(c1, Δ)
        dv, Δ = add_type_E(m, ens, c1)
        push!(step_list, Δ)
        apply_step!(c1, Δ)
        println("\n\nkinks1:")
        for i in c1.kinks
            println(i)
        end
        println("\n\nkinks2:\n\n\n")
        for i in apply_step(c2,step_list).kinks
            println(i)
        end
        dv, Δ = add_type_B(m, ens, c1)
        push!(step_list, Δ)
        apply_step!(c1, Δ)
    end
    println("\n\nkinks1:")
    for i in c1.kinks
        println(i)
    end
    println("\n\nkinks2:")
    for i in apply_step(c2,step_list).kinks
        println(i)
    end
    println("\n\n")
    #@test apply_step(c2,step_list).kinks == c1.kinks
    #@test apply_step(c2,step_list).occupations == c1.occupations
    #@test apply_step(c2,step_list) == c1
end
=#
