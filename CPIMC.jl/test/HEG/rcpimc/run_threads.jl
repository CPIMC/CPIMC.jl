using OnlineStats
using DelimitedFiles


include("../../../src/Configuration.jl")
include("../../../src/HEG/model.jl")
include("../../../src/HEG/RCPIMC/updates.jl")
include("../../../src/HEG/RCPIMC/estimators.jl")
include("../../../src/CPIMC.jl")


"""include("src/Configuration.jl")
include("src/HEG/model.jl")
include("src/CPIMC.jl")
include("src/HEG/RCPIMC/updates.jl")
include("src/HEG/RCPIMC/estimators.jl")"""
function main()
    # MC options
    NMC = 10^4
    cyc = 20
    NEquil = 10^4

    # system parameters
    θ = 1.0
    rs = 1.0
    S = get_sphere_with_same_spin(OrbitalHEG((0,0,0),1),dk=2)

    #4Particles
    #S = Set{Orbital{3}}([OrbitalHEG((0,0,0),1), OrbitalHEG((1,0,0),1), OrbitalHEG((0,1,0),1), OrbitalHEG((0,0,1),1)])

    println("#################################################")
    println("N: ", length(S))
    println("θ: ", θ)
    println("rs: ", rs)
    N = length(S)
    c = Configuration(S)

    e = Ensemble(rs, get_β_internal(θ,N), N) # get_β_internal only works for 3D
    updates = [move_particle, add_type_B, remove_type_B, change_type_B, shuffle_indices]#

    measurements = Dict(
      :Ekin => (Variance(), Ekin)
    , :W_off_diag => (Variance(), W_off_diag)
    , :W_diag => (Variance(), W_diag)
    , :Epot => (Variance(), Epot)
    , :K => (Variance(), K)
    , :E => (Variance(), E)
    , :occs => (Group([Variance() for i in 1:100]), occupations)
    )

    measurements_Mean = Dict(
      :Ekin => (Mean(), Ekin)
    , :W_off_diag => (Mean(), W_off_diag)
    , :W_diag => (Mean(), W_diag)
    , :Epot => (Mean(), Epot)
    , :K => (Mean(), K)
    , :E => (Mean(), E)
    , :occs => (Group([Mean() for i in 1:100]), occupations)
    )


    println("Start MC process ... ")
    Marcov_Chain_builders = Array{Task}(undef,Threads.nthreads())#die Anzahl threads ist inital die Anzahl Kerne
    Measurements_of_runs = Set{Dict{Symbol,Tuple{OnlineStat,Function}}}()
    for t in 1:Threads.nthreads()
        m = deepcopy(measurements_Mean)
        push!(Measurements_of_runs,m)
        Marcov_Chain_builders[t] = Threads.@spawn(runMC_multithreaded(NMC, cyc, NEquil, updates, m, e, c))
    end
    for mcb in Marcov_Chain_builders
        wait(mcb)
    end
    println(" finished.")
    println("measurements:")
    println("=============")
    #Avarage over the uncorrolated mean-values of the single runs
    for m in Measurements_of_runs
        for (key,(stat,obs)) in m
            if key == :occs
                fit!(first(measurements[key]), mean.(m[:occs][1].stats))
            else
                fit!(first(measurements[key]), mean(stat))
            end
        end
    end
    #print measurements
    for (k,(f,m)) in measurements
        if typeof(f) == Variance{Float64,Float64,EqualWeight}
            println(typeof(m).name.mt.name, "\t", mean(f), " +/- ", std(f)/sqrt(Threads.nthreads()-1))
            if  (k == :Epot) | (k == :E)
                println(typeof(m).name.mt.name,"_t_Ha", "\t", E_Ry(mean(f)-abs_E_mad(e.N, lambda(e.N,e.rs)),lambda(e.N,e.rs))/2, " +/- ", E_Ry(std(f),lambda(e.N,e.rs))/sqrt(Threads.nthreads()-1)/2)
            elseif (k == :Ekin)
                println(typeof(m).name.mt.name,"_Ha", "\t", E_Ry(mean(f),lambda(e.N,e.rs))/2, " +/- ", E_Ry(std(f),lambda(e.N,e.rs))/sqrt(Threads.nthreads()-1)/2)
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
    #open("test/HEG/rcpimc/out/occNums_$(N)_th$(replace(string(θ),"." => ""))_rs$(replace(string(rs),"." => ""))_Samples$((NMC*Threads.nthreads()/cyc)).dat", "w") do io
    open("CPIMC.jl/test/HEG/rcpimc/out/occNums_$(N)_th$(replace(string(θ),"." => ""))_rs$(replace(string(rs),"." => ""))_Samples$((NMC*Threads.nthreads()/cyc)).dat", "w") do io
           writedlm(io, zip(mean.(measurements[:occs][1].stats), std.(measurements[:occs][1].stats)/(NMC*Threads.nthreads()/cyc)))
    end
    #create resultsfile
    #open("test/HEG/rcpimc/out/results_$(N)_th$(replace(string(θ),"." => ""))_rs$(replace(string(rs),"." => ""))_Samples$((NMC*Threads.nthreads()/cyc)).dat", "w") do io
    open("CPIMC.jl/test/HEG/rcpimc/out/results_$(N)_th$(replace(string(θ),"." => ""))_rs$(replace(string(rs),"." => ""))_Samples$((NMC*Threads.nthreads()/cyc)).dat", "w") do io
        for (k,(f,m)) in measurements
            if typeof(f) == Variance{Float64,Float64,EqualWeight}
                write(io, string(typeof(m).name.mt.name, "\t", mean(f), " +/- ", std(f)/sqrt(Threads.nthreads()-1),"\n"))
                if  (k == :Epot) | (k == :E)
                    write(io, string(typeof(m).name.mt.name,"_t_Ha", "\t", E_Ry(mean(f)-abs_E_mad(e.N, lambda(e.N,e.rs)),lambda(e.N,e.rs))/2, " +/- ", (E_Ry(std(f),lambda(e.N,e.rs))/sqrt(Threads.nthreads()-1)/2),"\n"))
                elseif (k == :Ekin)
                    write(io,string(typeof(m).name.mt.name,"_Ha", "\t", E_Ry(mean(f),lambda(e.N,e.rs))/2, " +/- ", E_Ry(std(f),lambda(e.N,e.rs))/sqrt(Threads.nthreads()-1)/2,"\n"))
                end
            end
        end
    end
end

#Juno.@run(main())
main()
