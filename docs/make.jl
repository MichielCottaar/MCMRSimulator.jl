using MRSimulator
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(MRSimulator, :DocTestSetup, :(using MRSimulator); recursive=true)
bib = CitationBibliography(joinpath(@__DIR__, "references.bib"), sorting=:nyt)

makedocs(
    bib;
    modules=[MRSimulator],
    authors="Michiel Cottaar <Michiel.cottaar@ndcn.ox.ac.uk>",
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
        "Sequence" => "sequence.md",
        "Magnetic susceptibility" => "off_resonance.md",
        "API" => "api.md",
        "References" => "references.md",
    ],
)
