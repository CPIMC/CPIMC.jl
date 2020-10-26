using OnlineStats
using DelimitedFiles


include("../../../src/Configuration.jl")
include("../../../src/HEG/model.jl")
include("../../../src/HEG/RCPIMC/updates.jl")
include("../../../src/HEG/RCPIMC/estimators.jl")
include("../../../src/CPIMC.jl")


"""include("src/Configuration.jl")
include("src/CPIMC.jl")
include("src/HEG/model.jl")
include("src/HEG/RCPIMC/updates.jl")
include("src/HEG/RCPIMC/estimators.jl")"""

function main()
    # MC options
    NMC = 10^7
    cyc = 10
    NEquil = 10^6

    # system parameters
    theta = 0.5
    rs = 10

    #S = get_orbs_with_spin(get_sphere(Orbital((0,0,0),1),dk=2),1)

    #4Particles
    S = Set{Orbital{3}}([Orbital((0,0,0),1), Orbital((1,0,0),1), Orbital((0,1,0),1), Orbital((0,0,1),1)])

    println("Number of particles: ", length(S))
    println("theta: ", theta)
    println("rs: ", rs)
    N = length(S)
    println("N: ", N)
    c = Configuration(S)

    e = Ensemble(rs, get_beta_internal(theta,N), N) # get_beta_internal only works for 3D
    updates = [move_particle, Add_Type_B, remove_Type_B, change_type_B,shuffle_indixes]

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
    Marcov_Chain_builders = Array{Task}(undef,Threads.nthreads())#die Anzahl threads ist inital die Anzahl Kerne
    for t in 1:Threads.nthreads()
        Marcov_Chain_builders[t] = Threads.@spawn(runMC_multythreaded(NMC, cyc, NEquil, updates, measurements, e, c))
    end
    for mcb in Marcov_Chain_builders
        wait(mcb)
    end
    println(" finished.")
    println("measurements:")
    println("=============")

    for (k,(f,m)) in measurements
        if typeof(f) == Variance{Float64,Float64,EqualWeight}
            println(typeof(m).name.mt.name, "\t", mean(f), " +/- ", std(f))
            if  (k == :Epot) | (k == :E)
                println(typeof(m).name.mt.name,"_t_Ha", "\t", E_Ry(mean(f)-abs_E_Mad(e.N, lambda(e.N,e.rs)),lambda(e.N,e.rs))/2, " +/- ", E_Ry(std(f),lambda(e.N,e.rs))/2)
            elseif (k == :Ekin)
                println(typeof(m).name.mt.name,"_Ha", "\t", E_Ry(mean(f),lambda(e.N,e.rs))/2, " +/- ", E_Ry(std(f),lambda(e.N,e.rs))/2)
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
    #open("../out/occNums_N$(N)_th$(replace(string(theta),"." => ""))_rs$(replace(string(rs),"." => "")).dat", "w") do io
    #       writedlm(io, zip(mean.(measurements[:occs][1].stats), std.(measurements[:occs][1].stats)))
    #   end
end

#Juno.@run(main())
main()
