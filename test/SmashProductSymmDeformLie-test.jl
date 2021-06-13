lie = PD.lie
mod = PD.mod

@testset ExtendedTestSet "All PBWDeformations.SmashProductSymmDeformLie tests" begin
    @testset "smashProducSymmDeformLie constructor" begin
        @testset "A_2 with hw [1,1]" begin
            sp = PD.smashProductLie('A', 2, [1,1])
            deform = PD.smashProductSymmDeformLie('A', 2, [1,1])
            @test_broken deform.extraData.sp == sp.extraData
            @test deform.basis == sp.basis

            # Test that the module basis commutes and the other commutators come from the smash product
            @test deform.commTable == Dict(union(sp.commTable, [(mod(i),mod(j)) => [] for i in 1:sp.extraData.nV for j in 1:i-1]))
        end

        @testset "B_2 with hw [1,0]" begin
            sp = PD.smashProductLie('B', 2, [1,0])
            deform = PD.smashProductSymmDeformLie('B', 2, [1,0])
            @test_broken deform.extraData.sp == sp.extraData
            @test deform.basis == sp.basis

            # Test that the module basis commutes and the other commutators come from the smash product
            @test deform.commTable == Dict(union(sp.commTable, [(mod(i),mod(j)) => [] for i in 1:sp.extraData.nV for j in 1:i-1]))
        end
    end
    
end
