@testset "LinearIndependence.jl tests" begin
    @testset "is_linearly_independent(_with_relations)" begin
        #TODO

    end

    @testset "is_in_span" begin
        @test is_in_span(QQ, matrix(QQ, 2, 2, [1, 2, 3, 4]), [matrix(QQ, 2, 2, [1, 2, 3, 4])])

        @test !is_in_span(QQ, matrix(QQ, 2, 2, [1, 2, 3, 4]), [matrix(QQ, 2, 2, [1, 1, 1, 1])])

        @test is_in_span(
            QQ,
            matrix(QQ, 2, 2, [1, 2, 3, 4]),
            [matrix(QQ, 2, 2, [1, 2, 3, 4]), matrix(QQ, 2, 2, [0, 0, 0, 1])],
        )

        @test !is_in_span(
            QQ,
            matrix(QQ, 2, 2, [1, 2, 3, 4]),
            [matrix(QQ, 2, 2, [1, 2, 3, 6]), matrix(QQ, 2, 2, [1, 1, 1, 1])],
        )

        @test is_in_span(
            QQ,
            matrix(QQ, 2, 2, [1, 2, 3, 4]),
            [matrix(QQ, 2, 2, [1, 2, 3, 6]), matrix(QQ, 2, 2, [0, 0, 0, 1])],
        )

    end

    @testset "is_span_subset" begin
        @test is_span_subset(
            QQ,
            [matrix(QQ, 2, 2, [1, 2, 3, 4])],
            [matrix(QQ, 2, 2, [1, 2, 3, 6]), matrix(QQ, 2, 2, [0, 0, 0, 1])],
        )

        @test is_span_subset(
            QQ,
            [matrix(QQ, 2, 2, [1, 2, 3, 4]), matrix(QQ, 2, 2, [1, 2, 3, 5])],
            [matrix(QQ, 2, 2, [1, 2, 3, 6]), matrix(QQ, 2, 2, [0, 0, 0, 1])],
        )

        @test is_span_subset(
            QQ,
            [matrix(QQ, 2, 2, [1, 2, 3, 6]), matrix(QQ, 2, 2, [0, 0, 0, 1])],
            [matrix(QQ, 2, 2, [1, 2, 3, 4]), matrix(QQ, 2, 2, [1, 2, 3, 5])],
        )

        @test !is_span_subset(
            QQ,
            [matrix(QQ, 2, 2, [1, 2, 3, 6]), matrix(QQ, 2, 2, [0, 0, 0, 1])],
            [matrix(QQ, 2, 2, [1, 2, 3, 4]), matrix(QQ, 2, 2, [1, 1, 1, 1])],
        )

        @test !is_span_subset(
            QQ,
            [matrix(QQ, 2, 2, [1, 2, 3, 4]), matrix(QQ, 2, 2, [1, 1, 1, 1])],
            [matrix(QQ, 2, 2, [1, 2, 3, 6]), matrix(QQ, 2, 2, [0, 0, 0, 1])],
        )

        @test is_span_subset(
            QQ,
            [matrix(QQ, 2, 2, [1, 2, 3, 4]), matrix(QQ, 2, 2, [2, 4, 6, 8])],
            [matrix(QQ, 2, 2, [1, 2, 3, 4]), matrix(QQ, 2, 2, [1, 1, 1, 1])],
        )

        @test !is_span_subset(
            QQ,
            [matrix(QQ, 2, 2, [1, 2, 3, 4]), matrix(QQ, 2, 2, [1, 1, 1, 1])],
            [matrix(QQ, 2, 2, [1, 2, 3, 4]), matrix(QQ, 2, 2, [2, 4, 6, 8])],
        )

    end

    @testset "is_span_equal" begin
        @test !is_span_equal(
            QQ,
            [matrix(QQ, 2, 2, [1, 2, 3, 4]), matrix(QQ, 2, 2, [1, 1, 1, 1])],
            [matrix(QQ, 2, 2, [1, 2, 3, 6]), matrix(QQ, 2, 2, [0, 0, 0, 1])],
        )

        @test !is_span_equal(
            QQ,
            [matrix(QQ, 2, 2, [1, 2, 3, 4]), matrix(QQ, 2, 2, [1, 1, 1, 1])],
            [matrix(QQ, 2, 2, [1, 2, 3, 4]), matrix(QQ, 2, 2, [2, 4, 6, 8])],
        )

        @test is_span_equal(
            QQ,
            [matrix(QQ, 2, 2, [1, 2, 3, 4]), matrix(QQ, 2, 2, [1, 2, 3, 5])],
            [matrix(QQ, 2, 2, [1, 2, 3, 6]), matrix(QQ, 2, 2, [0, 0, 0, 1])],
        )

        @test is_span_equal(
            QQ,
            [matrix(QQ, 2, 2, [1, 2, 3, 4]), matrix(QQ, 2, 2, [1, 2, 3, 5]), matrix(QQ, 2, 2, [1, 2, 3, 6])],
            [matrix(QQ, 2, 2, [1, 2, 3, 6]), matrix(QQ, 2, 2, [0, 0, 0, 1])],
        )

    end


end
