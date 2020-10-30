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
    NMC = 10^6
    cyc = 100
    NEquil = 10^5

    # system parameters
    theta = 0.05
    rs = 1

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
            println(typeof(m).name.mt.name, "\t", mean(f), " +/- ", std(f)/sqrt(NMC*Threads.nthreads()/cyc))
            if  (k == :Epot) | (k == :E)
                println(typeof(m).name.mt.name,"_t_Ha", "\t", E_Ry(mean(f)-abs_E_Mad(e.N, lambda(e.N,e.rs)),lambda(e.N,e.rs))/2, " +/- ", E_Ry(std(f),lambda(e.N,e.rs))/sqrt(NMC*Threads.nthreads()/cyc)/2)
            elseif (k == :Ekin)
                println(typeof(m).name.mt.name,"_Ha", "\t", E_Ry(mean(f),lambda(e.N,e.rs))/2, " +/- ", E_Ry(std(f),lambda(e.N,e.rs))/sqrt(NMC*Threads.nthreads()/cyc)/2)
            end
        end
    end

    println("")

    #occupations funktionieren noch nicht fürs WW-System
    println("occupations:")
    println("============")
    println(mean.(measurements[:occs][1].stats))
    println(std.(measurements[:occs][1].stats))

    # create occnumsfile
    open("test/HEG/rcpimc/out/occNums_$(N)_th$(replace(string(theta),"." => ""))_rs$(replace(string(rs),"." => ""))_Samples$((NMC*Threads.nthreads()/cyc)).dat", "w") do io
           writedlm(io, zip(mean.(measurements[:occs][1].stats), std.(measurements[:occs][1].stats)/(NMC*Threads.nthreads()/cyc)))
    end
    #create resultsfile
    open("test/HEG/rcpimc/out/results_$(N)_th$(replace(string(theta),"." => ""))_rs$(replace(string(rs),"." => ""))_Samples$((NMC*Threads.nthreads()/cyc)).dat", "w") do io
        for (k,(f,m)) in measurements
            if typeof(f) == Variance{Float64,Float64,EqualWeight}
                write(io, string(typeof(m).name.mt.name, "\t", mean(f), " +/- ", std(f)/sqrt(NMC*Threads.nthreads()/cyc),"\n"))
                if  (k == :Epot) | (k == :E)
                    write(io, string(typeof(m).name.mt.name,"_t_Ha", "\t", E_Ry(mean(f)-abs_E_Mad(e.N, lambda(e.N,e.rs)),lambda(e.N,e.rs))/2, " +/- ", (E_Ry(std(f),lambda(e.N,e.rs))/sqrt(NMC*Threads.nthreads()/cyc)/2),"\n"))
                elseif (k == :Ekin)
                    write(io,string(typeof(m).name.mt.name,"_Ha", "\t", E_Ry(mean(f),lambda(e.N,e.rs))/2, " +/- ", E_Ry(std(f),lambda(e.N,e.rs))/sqrt(NMC*Threads.nthreads()/cyc)/2,"\n"))
                end
            end
        end
    end
end

#Juno.@run(main())
main()
