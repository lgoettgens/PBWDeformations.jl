using Documenter
using DocumenterCitations
using PBWDeformations

#! format: off

DocMeta.setdocmeta!(
    PBWDeformations,
    :DocTestSetup,
    :(using PBWDeformations; using PBWDeformations.Oscar);
    recursive = true,
)

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"); style=:alpha)

makedocs(
    bib,
    modules = [PBWDeformations],
    repo = "https://github.com/lgoettgens/PBWDeformations.jl/blob/{commit}{path}#{line}",
    sitename = "PBWDeformations.jl",
    checkdocs = :none, #:all, :exports
    strict = true,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://lgoettgens.github.io/PBWDeformations.jl/",
    ),
    pages = [
        "PBWDeformations.jl" => "index.md",
        "Smash products" => "smash_product_lie.md",
        "Smash product deformations" => "smash_product_deform_lie.md",
        "PBWDeformations" => "smash_product_pbwdeform_lie.md",
        "Arc diagrams" => "arc_diagrams.md",
        "Pseudographs" => "pseudographs.md",
        "Util functions" => "util.md",
        "References" => "references.md",
    ],
    # doctestfilters = [r"(Nemo\.)?QQFieldElem"],
)

deploydocs(
    repo = "github.com/lgoettgens/PBWDeformations.jl.git",
    deploy_config = Documenter.GitHubActions(),
    branch = "gh-pages",
)
