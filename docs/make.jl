using Documenter
using PBWDeformations

DocMeta.setdocmeta!(
    PBWDeformations,
    :DocTestSetup,
    :(using PBWDeformations);
    recursive = true,
)

makedocs(;
    modules = [PBWDeformations],
    repo = "https://gitlab.com/user/project/blob/{commit}{path}#{line}",
    sitename = "PBWDeformations.jl",
    checkdocs = :none,
    strict = true,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://johannesflake.gitlab.io/PBWDeformations.jl",
    ),
    pages = ["Home" => "index.md"],
)
