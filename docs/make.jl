using Documenter
using PBWDeformations

DocMeta.setdocmeta!(
    PBWDeformations,
    :DocTestSetup,
    :(using PBWDeformations; using Oscar);
    recursive = true,
)

makedocs(
    modules = [PBWDeformations],
    repo = "https://gitlab.com/johannesflake/pbwdeformations.jl/blob/{commit}{path}#{line}",
    sitename = "PBWDeformations.jl",
    checkdocs = :none,
    strict = true,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://johannesflake.gitlab.io/PBWDeformations.jl",
    ),
    pages = ["Home" => "index.md"],
    doctestfilters = [r"(Nemo\.)?fmpq"],
)

deploydocs(
    repo = "gitlab.com/johannesflake/pbwdeformations.jl.git",
    deploy_config = Documenter.GitLab(),
    branch = "pages",
    dirname = "public",
    devbranch = "master",
)
