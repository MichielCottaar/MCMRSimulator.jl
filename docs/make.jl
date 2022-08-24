using MRSimulator
using Documenter

DocMeta.setdocmeta!(MRSimulator, :DocTestSetup, :(using MRSimulator); recursive=true)

makedocs(;
    modules=[MRSimulator],
    authors="Michiel Cottaar <Mmichiel.cottaar@ndcn.ox.ac.uk>",
    repo="https://git.fmrib.ox.ac.uk/ndcn0236/MRSimulator.jl/blob/{commit}{path}#{line}",
    sitename="MRSimulator.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ndcn0236.gitlab.io/MRSimulator.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Geometry" => "obstructions.md",
        "Magnetic susceptibility" => "off_resonance.md",
        "API" => "api.md",
    ],
)
