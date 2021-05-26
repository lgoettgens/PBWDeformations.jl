using PBWDeformations
using Documenter

DocMeta.setdocmeta!(PBWDeformations, :DocTestSetup, :(using PBWDeformations); recursive=true)

makedocs(;
    modules=[PBWDeformations],
    authors="Johannes Flake <flake@art.rwth-aachen.de>, Lars Göttgens <lars.goettgens@rwth-aachen.de>, Phil Pützstück <phil.puetzstueck@rwth-aachen.de>",
    repo="https://gitlab.com/a-rt/PBWDeformations.jl/blob/{commit}{path}#{line}",
    sitename="PBWDeformations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://a-rt.gitlab.io/PBWDeformations.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    checkdocs=:exports,
    strict=true,
)
