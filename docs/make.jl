using MCMRSimulator
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(MCMRSimulator, :DocTestSetup, :(using MCMRSimulator); recursive=true)
bib = CitationBibliography(joinpath(@__DIR__, "references.bib"), sorting=:nyt)

makedocs(
    bib;
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
        "Geometry" => "obstructions.md",
        "Sequence" => "sequence.md",
        "Magnetic susceptibility" => "off_resonance.md",
        "API" => "api.md",
        "References" => "references.md",
    ],
    strict=[:example_block],
)
