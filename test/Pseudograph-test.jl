@testset ExtendedTestSet "All Pseudograph.jl tests" begin
    @testset "Pseudograph constructor" begin
        # negative weights
        for _ in 1:num_random_tests
            nloops = rand(0:10)
            nedges = rand(0:10)
            loops1 = rand(-2:20, nloops)
            loops2 = rand(-2:20, nloops)
            edges = rand(-2:20, nedges)

            if nloops > 0
                tmp = copy(loops1)
                tmp[rand(1:nloops)] = rand(-5:-1)
                @test_throws ArgumentError Pseudograph2(tmp, loops2, edges)

                tmp = copy(loops2)
                tmp[rand(1:nloops)] = rand(-5:-1)
                @test_throws ArgumentError Pseudograph2(loops1, tmp, edges)
            end

            if nedges > 0
                tmp = copy(edges)
                tmp[rand(1:nedges)] = rand(-5:-1)
                @test_throws ArgumentError Pseudograph2(loops1, loops2, tmp)
            end
        end

        # not regular
        for _ in 1:num_random_tests
            nloops1 = rand(0:10)
            nloops2 = rand(setdiff(0:10, nloops1))
            nedges = rand(0:10)
            loops1 = rand(0:20, nloops1)
            loops2 = rand(0:20, nloops2)
            edges = rand(0:20, nedges)

            @test_throws ArgumentError Pseudograph2(loops1, loops2, edges)
        end

        # correct
        @test Pseudograph2(Int[], Int[], Int[]) !== nothing
        @test Pseudograph2(Int[], Int[], Int[0]) !== nothing
        @test Pseudograph2(Int[], Int[], Int[1]) !== nothing
        @test Pseudograph2(Int[0], Int[0], Int[]) !== nothing
        @test Pseudograph2(Int[1], Int[0], Int[]) !== nothing
        for _ in 1:num_random_tests
            nloops = rand(0:10)
            nedges = rand(0:10)
            loops1 = rand(0:20, nloops)
            loops2 = rand(0:20, nloops)
            edges = rand(0:20, nedges)

            @test Pseudograph2(loops1, loops2, edges) !== nothing
        end
    end

    @testset "nedges and sum" begin
        for _ in 1:num_random_tests
            nloops = rand(0:10)
            ne = rand(0:10)
            loops1 = rand(0:20, nloops)
            loops2 = rand(0:20, nloops)
            edges = rand(0:20, ne)

            pg = Pseudograph2(loops1, loops2, edges)
            @test nedges(pg) == 2 * nloops + ne == length(loops1) + length(loops2) + length(edges)
            @test sum(pg) == sum(loops1) + sum(loops2) + sum(edges)
        end
    end

    @testset "all_pseudographs" begin
        @test length(all_pseudographs(0, 0)) == 1
        for s in 1:30
            @test length(all_pseudographs(0, s)) == 0
            @test length(all_pseudographs(1, s)) == 1
            @test length(all_pseudographs(2, s)) == div(3(s + 1) + 1, 2) # https://oeis.org/A007494
        end
    end

    @testset "all_pseudographs(upto_iso=true)" begin
        upto_iso = true
        @test length(all_pseudographs(0, 0; upto_iso)) == 1
        for s in 1:30
            @test length(all_pseudographs(0, s; upto_iso)) == 0
            @test length(all_pseudographs(1, s; upto_iso)) == 1
            @test length(all_pseudographs(2, s; upto_iso)) == 2 * div(s + 2, 2) # https://oeis.org/A052928
        end
    end

end
