pbw_arc_diagrams__so_extpowers_stdmod = PD.pbw_arc_diagrams__so_extpowers_stdmod

@testset ExtendedTestSet "All DeformationBases/*.jl tests" begin
    @testset "All DeformationBases/ArcDiagDeformBasis.jl tests" begin
        @testset "arcdiag_to_basiselem__so_extpowers_stdmod" begin
            sp, _ = smash_product_lie_so_extpowers_standard_module(QQ, 4, 2)
            @testset "not all specialisations are zero" begin
                diag = ArcDiagram(4, 2, [5, 3, 2, 6, 1, 4])
                dm = PD.arcdiag_to_basiselem__so_extpowers_stdmod(diag, 4, 2, 1, sp.alg(0), sp.basisL)
                @test !iszero(dm)
            end

            @testset "deformation is equivariant, d = $d" for d in 1:3
                for diag in pbw_arc_diagrams__so_extpowers_stdmod(2, d)
                    dm = PD.arcdiag_to_basiselem__so_extpowers_stdmod(diag, 4, 2, d, sp.alg(0), sp.basisL)
                    deform, _ = smash_product_deform_lie(sp, dm)
                    @test all(iszero, pbwdeform_eqs(deform, disabled=[:b, :c, :d]))
                end
            end
        end

    end

end
