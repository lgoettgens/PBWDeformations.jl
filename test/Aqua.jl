using Aqua

@testset "Aqua.jl" begin
    Aqua.test_all(
        PBWDeformations;
        ambiguities=false,          # recursive=false does not work
        unbound_args=true,
        undefined_exports=true,
        project_extras=true,
        stale_deps=true,
        deps_compat=true,
        project_toml_formatting=true,
        piracy=true,
    )
end
