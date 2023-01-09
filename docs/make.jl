using Documenter
using DocumenterCitations
using PBWDeformations

#! format: off

DocMeta.setdocmeta!(
    PBWDeformations,
    :DocTestSetup,
    :(using PBWDeformations; using Oscar);
    recursive = true,
)

bib = CitationBibliography("docs/references.bib", sorting = :nyt)

makedocs(
    bib,
    modules = [PBWDeformations],
    repo = "https://github.com/PBWDeformations/pbwdeformations.jl/blob/{commit}{path}#{line}",
    sitename = "PBWDeformations.jl",
    checkdocs = :none, #:all, :exports
    strict = true,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://pbwdeformations.github.io/pbwdeformations.jl/",
    ),
    pages = [
        "PBWDeformations.jl" => "index.md",
        "Smash products" => "smash_product_lie.md",
        "Smash product deformations" => "smash_product_deform_lie.md",
        "PBWDeformations" => "smash_product_pbwdeform_lie.md",
        "Arc diagrams" => "arc_diagrams.md",
        "Pseudographs" => "pseudographs.md",
        "Structure constants" => "structure_constants.md",
        "Util functions" => "util.md",
        "References" => "references.md",
    ],
    doctestfilters = [r"(Nemo\.)?fmpq"],
)

deploydocs(
    repo = "github.com/PBWDeformations/pbwdeformations.jl.git",
    deploy_config = Documenter.GitHubActions(),
    branch = "gh-pages",
)
