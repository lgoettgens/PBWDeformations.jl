@testset "liealgebra_sp" begin
    @testset "lie algebra; n = $n" for n in 1:6
        dimL = 2 * n^2 + n

        @test dimL == length(liealgebra_sp_basis(n))
        @test dimL == length(liealgebra_sp_symbols(n))
        @test (dimL, dimL) == size(liealgebra_sp_struct_const(n))
    end

    @testset "symmpowers std module; n = $n, e = $e" for n in 1:6, e in 1:min(5, n)
        dimL = length(liealgebra_sp_basis(n))
        dimV = binomial(2n + e - 1, e)
        @test dimV == length(liealgebra_sp_symmpowers_standard_module_symbols(n, e))
        @test (dimL, dimV) == size(liealgebra_sp_symmpowers_standard_module_struct_const(n, e))
    end

    @testset "extpowers std module; n = $n, e = $e" for n in 1:6, e in 1:min(5, n)
        dimL = length(liealgebra_sp_basis(n))
        dimV = binomial(2n, e)
        @test dimV == length(liealgebra_sp_extpowers_standard_module_symbols(n, e))
        @test (dimL, dimV) == size(liealgebra_sp_extpowers_standard_module_struct_const(n, e))
    end

end
