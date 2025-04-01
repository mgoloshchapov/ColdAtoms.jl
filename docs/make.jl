using Documenter, ColdAtoms

push!(LOAD_PATH,"../src/")

makedocs(
    modules=[ColdAtoms],
    sitename="ColdAtoms.jl",
    authors="M.Y. Goloshchapov",
    repo = Remotes.GitHub("mgoloshchapov", "ColdAtoms.jl"),
    pages = [
        "Home" => "index.md"
    ]
    )

deploydocs(
    repo = "https://github.com/mgoloshchapov/ColdAtoms.jl",
)