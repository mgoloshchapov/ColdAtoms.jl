using Documenter, ColdAtoms

makedocs(
    sitename="ColdAtoms.jl",
    authors="M.Y. Goloshchapov",
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API" => "library.md"
    ]
    )

deploydocs(
    repo = "https://github.com/mgoloshchapov/ColdAtoms.jl",
)