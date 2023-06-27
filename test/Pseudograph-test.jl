@testset "Pseudograph.jl tests" begin
    @testset "pseudograph constructor" begin
        # not regular
        for _ in 1:num_random_tests
            nloops1 = rand(0:10)
            nloops2 = rand(setdiff(0:10, nloops1))
            nedges = rand(0:10)
            loop1_weights = rand(0:20, nloops1)
            loop2_weights = rand(0:20, nloops2)
            edge_weights = rand(0:20, nedges)
            edges = MSet(
                [
                    [MSet([1, 1]) => k for k in loop1_weights]
                    [MSet([2, 2]) => k for k in loop2_weights]
                    [MSet([1, 2]) => k for k in edge_weights]
                ],
            )
            @test_throws ArgumentError PseudographLabelled(2, edges; regular_degree=2 * nloops1 + nedges)
        end

        # correct
        @test PseudographLabelled(2, Pair{MSet{Int}, Int}[]) !== nothing
        @test PseudographLabelled(2, [MSet([1, 2]) => 0]) !== nothing
        @test PseudographLabelled(2, [MSet([1, 2]) => 0]; regular_degree=1) !== nothing
        @test PseudographLabelled(2, [MSet([1, 1]) => 1, MSet([2, 2]) => 0]) !== nothing
        @test PseudographLabelled(2, [MSet([1, 1]) => 1, MSet([2, 2]) => 0]; regular_degree=2) !== nothing
        for _ in 1:num_random_tests
            numloops = rand(0:10)
            numedges = rand(0:10)
            loop1_weights = rand(0:20, numloops)
            loop2_weights = rand(0:20, numloops)
            edge_weights = rand(0:20, numedges)
            edges = MSet(
                [
                    [MSet([1, 1]) => k for k in loop1_weights]
                    [MSet([2, 2]) => k for k in loop2_weights]
                    [MSet([1, 2]) => k for k in edge_weights]
                ],
            )
            @test PseudographLabelled(2, edges) !== nothing
        end
    end

    @testset "nedges and sum" begin
        for _ in 1:num_random_tests
            numloops = rand(0:10)
            numedges = rand(0:10)
            loop1_weights = rand(0:20, numloops)
            loop2_weights = rand(0:20, numloops)
            edge_weights = rand(0:20, numedges)
            edges = MSet(
                [
                    [MSet([1, 1]) => k for k in loop1_weights]
                    [MSet([2, 2]) => k for k in loop2_weights]
                    [MSet([1, 2]) => k for k in edge_weights]
                ],
            )
            pg = PseudographLabelled(2, edges)
            @test nedges(pg) ==
                  2 * numloops + numedges ==
                  length(loop1_weights) + length(loop2_weights) + length(edge_weights)
            @test sum(pg) == sum(loop1_weights) + sum(loop2_weights) + sum(edge_weights)
        end
    end

    @testset "all_pseudographs" begin
        @test length(all_pseudographs(2, 0, 0)) == 1
        for s in 1:30
            @test length(all_pseudographs(2, 0, s)) == 0
            @test length(all_pseudographs(2, 1, s)) == 1
            @test length(all_pseudographs(2, 2, s)) == div(3(s + 1) + 1, 2) # https://oeis.org/A007494
        end
    end

    @testset "all_pseudographs(upto_iso=true)" begin
        upto_iso = true
        @test length(all_pseudographs(2, 0, 0; upto_iso)) == 1
        for s in 1:30
            @test length(all_pseudographs(2, 0, s; upto_iso)) == 0
            @test length(all_pseudographs(2, 1, s; upto_iso)) == 1
            @test length(all_pseudographs(2, 2, s; upto_iso)) == 2 * div(s + 2, 2) # https://oeis.org/A052928
        end
    end

end
