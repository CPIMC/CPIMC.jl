#TODO Averaging over Differnet Runs with same parameters
using DataFrames
using CSV

Filenames = filter(x -> occursin(r"results_N.*", x), readdir("test/UEG/Full_CPIMC/out"))

Collected_Results = DataFrame()
for Filename in Filenames


    local df = DataFrame(CSV.File(string("test/UEG/Full_CPIMC/out/", Filename)))

    row_df = DataFrame(
    N = parse(Int,match(r"\d\d*", match(r"N\d\d*",Filename).match).match),
    θ = parse(Float64,match(r"\d*\.\d*", match(r"th\d*\.\d*",Filename).match).match),
    rs = parse(Float64,match(r"\d*\.\d*", match(r"rs\d*\.\d*",Filename).match).match),
    steps = parse(Float64,match(r"\d*\.\d*", match(r"steps\d*\.\d*",Filename).match).match),
    sign = df[!,:sign][1],
    sign_err = df[!,:sign_err][1],
    Ekin = df[!,:K_fermion][1],
    Ekin_err = df[!,:Ekin_err][1],
    W_diag = df[!,:K][1],
    W_diag_err = df[!,:K_err][1],
    K = df[!,:K][1],
    K_err = df[!,:K_err][1],
    K_fermion = df[!,:K_fermion][1],
    K_fermion_err = df[!,:K_fermion_err][1],
    W = df[!,:W][1],
    ΔW = df[!,:ΔW][1],
    E = df[!,:E][1],
    ΔE = df[!,:ΔE][1],
    Wt_Ry = df[!,:Wt_Ry][1],
    ΔWt_Ry = df[!,:ΔWt_Ry][1],
    T_Ry = df[!,:T_Ry][1],
    ΔT_Ry = df[!,:ΔT_Ry][1],
    E_Ry = df[!,:E_Ry][1],
    ΔE_Ry = df[!,:ΔE_Ry][1]
    )
    append!(Collected_Results,row_df)
end

CSV.write("test/UEG/Full_CPIMC/out/Collected_Results.csv",Collected_Results)
