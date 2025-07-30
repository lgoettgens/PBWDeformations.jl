@testset "Misc.jl tests" begin
    @testset "is_prefix_equal" begin
        is_prefix_equal = PBWDeformations.is_prefix_equal

        @test is_prefix_equal([1, 2, 3], [1, 2, 3])
        @test !is_prefix_equal([1, 2, 3], [1, 2, 4])
        @test is_prefix_equal([1, 2], [1, 2, 3])
        @test is_prefix_equal([1, 2, 3], [1, 2])
        @test !is_prefix_equal([1, 2], [1, 3])
    end

    @testset "symmetrize" begin
        R, (x,y,z, t) = free_associative_algebra(QQ, [:x,:y,:z, :t])

        @test symmetrize(x) == x

        @test symmetrize(x*y) == QQ(1//2) * (x*y + y*x)
        @test symmetrize(y^2) == y^2

        @test symmetrize(x*y*z) == QQ(1//6) * (x*y*z + x*z*y + y*x*z + y*z*x + z*x*y + z*y*x)
        @test allequal([symmetrize(x*y*z), symmetrize(x*z*y), symmetrize(y*x*z), symmetrize(y*z*x), symmetrize(z*x*y), symmetrize(z*y*x)])
        @test symmetrize(x^2*y) == QQ(2//6) * (x^2*y + x*y*x + y*x^2)
        @test allequal([symmetrize(x^2*y), symmetrize(x*y*x), symmetrize(y*x^2)])
        @test symmetrize(z^3) == z^3

        @test symmetrize(x*y*z*t) == QQ(1//24) * (
            x*y*z*t + x*y*t*z + x*z*y*t + x*z*t*y + x*t*y*z + x*t*z*y +
            y*x*z*t + y*x*t*z + y*z*x*t + y*z*t*x + y*t*x*z + y*t*z*x +
            z*x*y*t + z*x*t*y + z*y*x*t + z*y*t*x + z*t*x*y + z*t*y*x +
            t*x*y*z + t*x*z*y + t*y*x*z + t*y*z*x + t*z*x*y + t*z*y*x
        )
        @test symmetrize(x^2*y*z) == QQ(2//24) * (
            x^2*y*z + x^2*z*y + x*y*x*z + x*y*z*x + x*z*x*y + x*z*y*x + y*x^2*z + y*x*z*x + y*z*x^2 + z*x^2*y + z*x*y*x + z*y*x^2
        )
        @test allequal([
            symmetrize(x^2 * y * z),
            symmetrize(x^2 * z * y),
            symmetrize(x * y * x * z),
            symmetrize(x * y * z * x),
            symmetrize(x * z * x * y),
            symmetrize(x * z * y * x),
            symmetrize(y * x^2 * z),
            symmetrize(y * x * z * x),
            symmetrize(y * z * x^2),
            symmetrize(z * x^2 * y),
            symmetrize(z * x * y * x),
            symmetrize(z * y * x^2),
        ])
        @test symmetrize(x^3*y) == QQ(6//24) * (x^3*y + x^2*y*x + x*y*x^2 + y*x^3)
        @test allequal([symmetrize(x^3*y), symmetrize(x^2*y*x), symmetrize(x*y*x^2), symmetrize(y*x^3)])
        @test symmetrize(x^2*y^2) == QQ(2*2//24) * (x*x*y*y + x*y*x*y + x*y*y*x + y*x*x*y + y*x*y*x + y*y*x*x)
        @test allequal([symmetrize(x*x*y*y), symmetrize(x*y*x*y), symmetrize(x*y*y*x), symmetrize(y*x*x*y), symmetrize(y*x*y*x), symmetrize(y*y*x*x)])
        @test symmetrize(t^4) == t^4

        # test linearity
        @test symmetrize(2*y^4*z*x - z*x*z^2*y*x + 9*x^2*z*x^2 + QQ(5//3)*y*z*x^2*z - 4*x*y^2 + z*x^2 + QQ(9//5)*z*x) ==
            2*symmetrize(y^4*z*x) - symmetrize(z*x*z^2*y*x) + 9*symmetrize(x^2*z*x^2) + QQ(5//3)*symmetrize(y*z*x^2*z) - 4*symmetrize(x*y^2) + symmetrize(z*x^2) + QQ(9//5)*symmetrize(z*x)
    end
end
