using Documenter, CPIMC

makedocs(
    modules = [CPIMC],
    sitename = "cpimc2020",
    authors = "Kai Hunger, Arif Yilmaz, Paul Hamann",
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "cpimc.md",
            "estimators.md",
            "planewaves.md",
            "ueg.md"
        ]])

deploydocs(repo="https://gitlab.physik.uni-kiel.de/hunger/cpimc2020")



