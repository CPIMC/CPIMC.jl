using Documenter, CPIMC

makedocs(
    modules = [CPIMC],
    sitename = "cpimc2020",
    authors = "Kai Hunger, Arif Yilmaz, Paul Hamann",
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "cpimc.md",
            "configuration.md",
            "estimators.md",
            "planewaves.md"
        ]])

deploydocs(repo="https://gitlab.physik.uni-kiel.de/hunger/cpimc2020")



