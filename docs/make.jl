using MCMRSimulator
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(MCMRSimulator, :DocTestSetup, :(using MCMRSimulator); recursive=true)
bib = CitationBibliography(joinpath(@__DIR__, "references.bib"), style=:authoryear)

makedocs(;
    modules=[MCMRSimulator],
    authors="Michiel Cottaar <Michiel.cottaar@ndcn.ox.ac.uk>",
    repo="https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/blob/{commit}{path}#{line}",
    sitename="MCMRSimulator.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial (Julia)" => "tutorial_julia.md",
        "Tutorial (CLI)" => "tutorial_cli.md",
        "Geometry" => "geometry.md",
        "Sequence" => "sequence.md",
        "MRI/collision properties" => "properties.md",
        "API" => "api.md",
        "References" => "references.md",
    ],
    warnonly=Documenter.except(:example_block),
    plugins=[bib],
)

deploydocs(repo="git.fmrib.ox.ac.uk:ndcn0236/mcmrsimulator.jl.git", branch="pages", devbranch="main")