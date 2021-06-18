lie = PD.lie
mod = PD.mod

@testset ExtendedTestSet "All PBWDeformations.SmashProductDeformLie tests" begin
    @testset "smashProducDeformLie coincides with smashProductSymmDeformLie on symmetric kappa" begin
        @testset "B_2 with hw [1,0]" begin
            sp = PD.smashProductLie('B', 2, [1,0])
            nV = sp.extraData.nV

            kappa = fill(PD.AlgebraElement(), nV, nV)

            deform1 = PD.smashProductSymmDeformLie(sp)
            deform2 = PD.smashProductDeformLie(sp, kappa)

            @test deform1 == deform2
        end
    end
    
end
