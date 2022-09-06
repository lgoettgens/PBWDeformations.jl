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
    repo = "https://gitlab.com/johannesflake/pbwdeformations.jl/blob/{commit}{path}#{line}",
    sitename = "PBWDeformations.jl",
    checkdocs = :none,
    strict = true,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://johannesflake.gitlab.io/PBWDeformations.jl",
    ),
    pages = [
        "PBWDeformations.jl" => "index.md",
        "Smash products" => "smash_product_lie.md",
        "Smash product deformations" => "smash_product_deform_lie.md",
        "Structure constants" => "structure_constants.md",
        "Util functions" => "util.md",
        "References" => "references.md",
    ],
    doctestfilters = [r"(Nemo\.)?fmpq"],
)

deploydocs(
    repo = "gitlab.com/johannesflake/pbwdeformations.jl.git",
    deploy_config = Documenter.GitLab(),
    branch = "pages",
    dirname = "public",
    devbranch = "master",
)
