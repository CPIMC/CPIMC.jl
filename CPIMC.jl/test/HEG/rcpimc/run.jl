using OnlineStats
using DelimitedFiles


include("../../../src/Configuration.jl")
include("../../../src/HEG/model.jl")
include("../../../src/HEG/RCPIMC/updates.jl")
include("../../../src/HEG/RCPIMC/estimators.jl")
include("../../../src/CPIMC.jl")


"""include("CPIMC.jl/src/Configuration.jl")
include("CPIMC.jl/src/HEG/model.jl")
include("CPIMC.jl/src/CPIMC.jl")
include("CPIMC.jl/src/HEG/RCPIMC/updates.jl")
include("CPIMC.jl/src/HEG/RCPIMC/estimators.jl")"""

function main()
    # MC options
    NMC = 10^5###############################################################
    cyc = 20
    NEquil = 10^5
    #auffälligerBalken um nicht übersehbaren unterschied im vergleich zu run_threads herzustellen
    """#####################################################################
    ########################################################################
    ####################################################################"""
    # system parameters
    θ = 0.5
    rs = 10

    S = get_sphere_with_same_spin(OrbitalHEG((0,0,0),1),dk=1)

    #4Particles
    #S = Set{OrbitalHEG{3}}([OrbitalHEG((0,0,0),1), OrbitalHEG((1,0,0),1), OrbitalHEG((0,1,0),1), OrbitalHEG((0,0,1),1)])

    println("Number of particles: ", length(S))
    println("θ: ", θ)
    println("rs: ", rs)
    N = length(S)
    println("N: ", N)
    c = Configuration(S)

    e = Ensemble(rs, get_β_internal(θ,N), N) # get_β_internal only works for 3D
    updates = [move_particle, add_type_B, remove_type_B, change_type_B,shuffle_indices]

    measurements = Dict(
      :Ekin => (Variance(), Ekin)
    , :W_off_diag => (Variance(), W_off_diag)
    , :W_diag => (Variance(), W_diag)
    , :Epot => (Variance(), Epot)
    , :K => (Variance(), K)
    , :E => (Variance(), E)
    , :occs => (Group([Variance() for i in 1:100]), occupations)
    )

    println("Start MC process ... ")
    runMC(NMC, cyc, NEquil, updates, measurements, e, c)
    println(" finished.")
    println("measurements:")
    println("=============")

    for (k,(f,m)) in measurements
        if typeof(f) == Variance{Float64,Float64,EqualWeight}
            println(typeof(m).name.mt.name, "\t", mean(f), " +/- ", std(f)/sqrt(NMC/cyc))
            if  (k == :Epot) | (k == :E)
                println(typeof(m).name.mt.name,"_t_Ha", "\t", E_Ry(mean(f)-abs_E_mad(e.N, lambda(e.N,e.rs)),lambda(e.N,e.rs))/2, " +/- ", E_Ry(std(f),lambda(e.N,e.rs))/sqrt(NMC/cyc)/2)
            elseif (k == :Ekin)
                println(typeof(m).name.mt.name,"_Ha", "\t", E_Ry(mean(f),lambda(e.N,e.rs))/2, " +/- ", E_Ry(std(f),lambda(e.N,e.rs))/sqrt(NMC/cyc)/2)
            end
        end
    end

    println("")

    #occupations funktionieren noch nicht fürs WW-System
    println("occupations:")
    println("============")
    println(mean.(measurements[:occs][1].stats))
    println(std.(measurements[:occs][1].stats))

    # Print to results file
    #open("../out/occNums_N$(N)_th$(replace(string(θ),"." => ""))_rs$(replace(string(rs),"." => "")).dat", "w") do io
    #       writedlm(io, zip(mean.(measurements[:occs][1].stats), std.(measurements[:occs][1].stats)))
    #   end
end

#Juno.@run(main())
main()
