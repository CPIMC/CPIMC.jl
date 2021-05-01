using Documenter, CPIMC

makedocs(
    modules = [CPIMC],
    sitename = "CPIMC.jl",
    authors = "Kai Hunger, Arif Yilmaz, Paul Hamann",
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "cpimc.md",
            "estimators.md",
            "planewaves.md",
            "ueg.md"
        ]])

deploydocs(
    repo = "github.com/CPIMC/CPIMC.jl.git",
)

