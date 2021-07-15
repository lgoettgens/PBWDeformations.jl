@testset ExtendedTestSet "All PBWDeformations.SmashProductDeformLie tests" begin
    @testset "smashProducDeformLie coincides with smashProductSymmDeformLie on symmetric kappa" begin
        @testset "B_2 with hw [1,0]" begin
            sp = PD.smashProductLie('B', 2, [1,0])
            nV = sp.extraData.nV

            kappa = fill(algebraElement(), nV, nV)

            deform1 = PD.smashProductSymmDeformLie(sp)
            deform2 = PD.smashProductDeformLie(sp, kappa)

            @test deform1 == deform2
        end
    end

    @testset "smashProducSymmDeformLie constructor" begin
        @testset "$(dynkin)_$n with hw $lambda" for (dynkin, n, lambda) in [('A', 2, [1,1]), ('B', 2, [1,0])]
            sp = PD.smashProductLie(dynkin, n, lambda)
            deform1 = PD.smashProductSymmDeformLie(sp)
            deform2 = PD.smashProductSymmDeformLie(dynkin, n, lambda)

            @test deform1.extraData.symmetric == deform2.extraData.symmetric == true
            @test typeof(deform1.extraData.sp) == typeof(deform2.extraData.sp) == typeof(sp.extraData)
            @test deform1.extraData.sp == deform2.extraData.sp == sp.extraData
            @test deform1.basis == deform2.basis == sp.basis

            # Test that the module basis commutes and the other commutators come from the smash product
            @test deform1.relTable == deform2.relTable == Dict(union(sp.relTable, [(mod(i), mod(j)) => [(1, [mod(j), mod(i)])] for i in 1:sp.extraData.nV for j in 1:i-1]))

            @test sprint(show, deform1) == sprint(show, deform2)

            showOutput = @test_nowarn sprint(show, deform1)
            @test occursin("smash product", lowercase(showOutput))
            @test occursin("lie algebra", lowercase(showOutput))
            @test occursin(string(dynkin), showOutput)
            @test occursin(string(lambda), showOutput)
            @test occursin("symmetric deformation", lowercase(showOutput))
        end

    end

    @testset "smashProductDeformLie assertions" begin
        @testset "assert correct dimensions of kappa" begin
            sp = PD.smashProductLie('B', 2, [1,0])
            nV = sp.extraData.nV

            kappa = fill(algebraElement(), nV+1, nV)
            @test_throws AssertionError("size of kappa matches module dimension") PD.smashProductDeformLie(sp, kappa)

            kappa = fill(algebraElement(), nV, nV+1)
            @test_throws AssertionError("size of kappa matches module dimension") PD.smashProductDeformLie(sp, kappa)

            kappa = fill(algebraElement(), nV+1, nV+1)
            @test_throws AssertionError("size of kappa matches module dimension") PD.smashProductDeformLie(sp, kappa)
        end

        @testset "assert entries of kappa contained in Hopf algebra of smash product" begin
            sp = PD.smashProductLie('B', 2, [1,0])
            nV = sp.extraData.nV

            # basis of sp consists of (:mod, 1) to (:mod 5) and (:lie, 1) to (:lie, 10)

            kappa = fill(algebraElement(), nV, nV)
            kappa[1,2] = algebraElement(mod(1))
            kappa[2,1] = -kappa[1,2]
            @test_throws AssertionError("kappa only takes values in Hopf algebra") PD.smashProductDeformLie(sp, kappa)

            kappa = fill(algebraElement(), nV, nV)
            kappa[1,2] = algebraElement(lie(0))
            kappa[2,1] = -kappa[1,2]
            @test_throws AssertionError("kappa only takes values in Hopf algebra") PD.smashProductDeformLie(sp, kappa)

            kappa = fill(algebraElement(), nV, nV)
            kappa[1,2] = algebraElement(lie(11))
            kappa[2,1] = -kappa[1,2]
            @test_throws AssertionError("kappa only takes values in Hopf algebra") PD.smashProductDeformLie(sp, kappa)
        end

        @testset "assert kappa is skew symmetric" begin
            sp = PD.smashProductLie('B', 2, [1,0])
            nV = sp.extraData.nV

            kappa = fill(PD.AlgebraElement(), nV, nV)
            kappa[1,1] = algebraElement(lie(1))
            @test_throws AssertionError("kappa is skew-symmetric") PD.smashProductDeformLie(sp, kappa)

            kappa = fill(PD.AlgebraElement(), nV, nV)
            kappa[1,2] = algebraElement(lie(1))
            @test_throws AssertionError("kappa is skew-symmetric") PD.smashProductDeformLie(sp, kappa)
        end

    end

end
