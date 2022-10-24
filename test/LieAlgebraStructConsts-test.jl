liealgebra_so_basis = PD.liealgebra_so_basis
liealgebra_so_symbols = PD.liealgebra_so_symbols
liealgebra_so_struct_const = PD.liealgebra_so_struct_const
liealgebra_so_standard_module_basis = PD.liealgebra_so_standard_module_basis
liealgebra_so_fundamental_module_symbols = PD.liealgebra_so_fundamental_module_symbols
liealgebra_so_fundamental_module_struct_const = PD.liealgebra_so_fundamental_module_struct_const
liealgebra_so_symmpowers_standard_module_symbols = PD.liealgebra_so_symmpowers_standard_module_symbols
liealgebra_so_symmpowers_standard_module_struct_const = PD.liealgebra_so_symmpowers_standard_module_struct_const
liealgebra_so_extpowers_standard_module_symbols = PD.liealgebra_so_extpowers_standard_module_symbols
liealgebra_so_extpowers_standard_module_struct_const = PD.liealgebra_so_extpowers_standard_module_struct_const
liealgebra_sp_basis = PD.liealgebra_sp_basis
liealgebra_sp_symbols = PD.liealgebra_sp_symbols
liealgebra_sp_struct_const = PD.liealgebra_sp_struct_const
liealgebra_sp_standard_module_basis = PD.liealgebra_sp_standard_module_basis
liealgebra_sp_symmpowers_standard_module_symbols = PD.liealgebra_sp_symmpowers_standard_module_symbols
liealgebra_sp_symmpowers_standard_module_struct_const = PD.liealgebra_sp_symmpowers_standard_module_struct_const
liealgebra_sp_extpowers_standard_module_symbols = PD.liealgebra_sp_extpowers_standard_module_symbols
liealgebra_sp_extpowers_standard_module_struct_const = PD.liealgebra_sp_extpowers_standard_module_struct_const

@testset ExtendedTestSet "All LieAlgebraStructConsts.jl tests" begin
    @testset "liealgebra_so" begin
        @testset "lie algebra; n = $n" for n in 1:7
            if n % 2 == 1
                l = div(n - 1, 2)
                dimL = 2 * l^2 + l
            else
                l = div(n, 2)
                dimL = 2 * l^2 - l
            end

            @test dimL == length(liealgebra_so_basis(n))
            @test dimL == length(liealgebra_so_symbols(n))
            @test (dimL, dimL) == size(liealgebra_so_struct_const(n))
        end

        @testset "symmpowers std module; n = $n, e = $e" for n in 1:7, e in 1:min(5, n)
            dimL = length(liealgebra_so_basis(n))
            dimV = binomial(n + e - 1, e)
            @test dimV == length(liealgebra_so_symmpowers_standard_module_symbols(n, e))
            @test (dimL, dimV) == size(liealgebra_so_symmpowers_standard_module_struct_const(n, e))
        end

        @testset "extpowers std module; n = $n, e = $e" for n in 1:7, e in 1:min(5, n)
            dimL = length(liealgebra_so_basis(n))
            dimV = binomial(n, e)
            @test dimV == length(liealgebra_so_extpowers_standard_module_symbols(n, e))
            @test (dimL, dimV) == size(liealgebra_so_extpowers_standard_module_struct_const(n, e))
        end

    end

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

end
