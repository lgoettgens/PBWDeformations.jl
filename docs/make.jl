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

warnonly = false

makedocs(;
    modules = [PBWDeformations],
    sitename = "PBWDeformations.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
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
    clean=true,
    warnonly=warnonly,
    treat_markdown_warnings_as_error=!warnonly,
    checkdocs=:none,
    plugins = [bib],
)

deploydocs(
    repo = "github.com/lgoettgens/PBWDeformations.jl.git",
    deploy_config = Documenter.GitHubActions(),
    branch = "gh-pages",
)
