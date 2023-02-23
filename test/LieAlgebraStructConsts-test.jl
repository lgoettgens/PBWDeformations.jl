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

        @testset "correctness regression" begin
            @test repr("text/plain", PD.liealgebra_so_struct_const(3)) ==
                  "3×3 Matrix{Vector{Tuple{Int64, Int64}}}:\n []         [(-1, 3)]  [(1, 2)]\n [(1, 3)]   []         [(-1, 1)]\n [(-1, 2)]  [(1, 1)]   []"
            @test repr("text/plain", PD.liealgebra_so_struct_const(4)) ==
                  "6×6 Matrix{Vector{Tuple{Int64, Int64}}}:\n []         [(-1, 4)]  [(-1, 5)]  [(1, 2)]   [(1, 3)]   []\n [(1, 4)]   []         [(-1, 6)]  [(-1, 1)]  []         [(1, 3)]\n [(1, 5)]   [(1, 6)]   []         []         [(-1, 1)]  [(-1, 2)]\n [(-1, 2)]  [(1, 1)]   []         []         [(-1, 6)]  [(1, 5)]\n [(-1, 3)]  []         [(1, 1)]   [(1, 6)]   []         [(-1, 4)]\n []         [(-1, 3)]  [(1, 2)]   [(-1, 5)]  [(1, 4)]   []"
            @test repr("text/plain", PD.liealgebra_so_symmpowers_standard_module_struct_const(3, 1)) ==
                  "3×3 Matrix{Vector{Tuple{Int64, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []\n [(-1, 3)]  []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]"
            @test repr("text/plain", PD.liealgebra_so_symmpowers_standard_module_struct_const(3, 2)) ==
                  "3×6 Matrix{Vector{Tuple{Int64, Int64}}}:\n [(-1, 2), (-1, 2)]  [(-1, 4), (1, 1)]  [(-1, 5)]          [(1, 2), (1, 2)]    [(1, 3)]           []\n [(-1, 3), (-1, 3)]  [(-1, 5)]          [(-1, 6), (1, 1)]  []                  [(1, 2)]           [(1, 3), (1, 3)]\n []                  [(-1, 3)]          [(1, 2)]           [(-1, 5), (-1, 5)]  [(-1, 6), (1, 4)]  [(1, 5), (1, 5)]"
            @test repr("text/plain", PD.liealgebra_so_symmpowers_standard_module_struct_const(3, 3)) ==
                  "3×10 Matrix{Vector{Tuple{Int64, Int64}}}:\n [(-1, 2), (-1, 2), (-1, 2)]  [(-1, 4), (-1, 4), (1, 1)]  [(-1, 5), (-1, 5)]          [(-1, 7), (1, 2), (1, 2)]  [(-1, 8), (1, 3)]  [(-1, 9)]                   [(1, 4), (1, 4), (1, 4)]     [(1, 5), (1, 5)]            [(1, 6)]                    []\n [(-1, 3), (-1, 3), (-1, 3)]  [(-1, 5), (-1, 5)]          [(-1, 6), (-1, 6), (1, 1)]  [(-1, 8)]                  [(-1, 9), (1, 2)]  [(-1, 10), (1, 3), (1, 3)]  []                           [(1, 4)]                    [(1, 5), (1, 5)]            [(1, 6), (1, 6), (1, 6)]\n []                           [(-1, 3)]                   [(1, 2)]                    [(-1, 5), (-1, 5)]         [(-1, 6), (1, 4)]  [(1, 5), (1, 5)]            [(-1, 8), (-1, 8), (-1, 8)]  [(-1, 9), (-1, 9), (1, 7)]  [(-1, 10), (1, 8), (1, 8)]  [(1, 9), (1, 9), (1, 9)]"
            @test repr("text/plain", PD.liealgebra_so_symmpowers_standard_module_struct_const(4, 1)) ==
                  "6×4 Matrix{Vector{Tuple{Int64, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []         []\n [(-1, 3)]  []         [(1, 1)]   []\n [(-1, 4)]  []         []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]   []\n []         [(-1, 4)]  []         [(1, 2)]\n []         []         [(-1, 4)]  [(1, 3)]"
            @test repr("text/plain", PD.liealgebra_so_extpowers_standard_module_struct_const(3, 1)) ==
                  "3×3 Matrix{Vector{Tuple{Int64, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []\n [(-1, 3)]  []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]"
            @test repr("text/plain", PD.liealgebra_so_extpowers_standard_module_struct_const(3, 2)) ==
                  "3×3 Matrix{Vector{Tuple{Int64, Int64}}}:\n []         [(-1, 3)]  [(1, 2)]\n [(1, 3)]   []         [(-1, 1)]\n [(-1, 2)]  [(1, 1)]   []"
            @test repr("text/plain", PD.liealgebra_so_extpowers_standard_module_struct_const(4, 1)) ==
                  "6×4 Matrix{Vector{Tuple{Int64, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []         []\n [(-1, 3)]  []         [(1, 1)]   []\n [(-1, 4)]  []         []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]   []\n []         [(-1, 4)]  []         [(1, 2)]\n []         []         [(-1, 4)]  [(1, 3)]"
            @test repr("text/plain", PD.liealgebra_so_extpowers_standard_module_struct_const(4, 2)) ==
                  "6×6 Matrix{Vector{Tuple{Int64, Int64}}}:\n []         [(-1, 4)]  [(-1, 5)]  [(1, 2)]   [(1, 3)]   []\n [(1, 4)]   []         [(-1, 6)]  [(-1, 1)]  []         [(1, 3)]\n [(1, 5)]   [(1, 6)]   []         []         [(-1, 1)]  [(-1, 2)]\n [(-1, 2)]  [(1, 1)]   []         []         [(-1, 6)]  [(1, 5)]\n [(-1, 3)]  []         [(1, 1)]   [(1, 6)]   []         [(-1, 4)]\n []         [(-1, 3)]  [(1, 2)]   [(-1, 5)]  [(1, 4)]   []"
            @test repr("text/plain", PD.liealgebra_so_extpowers_standard_module_struct_const(4, 3)) ==
                  "6×4 Matrix{Vector{Tuple{Int64, Int64}}}:\n []         []         [(-1, 4)]  [(1, 3)]\n []         [(1, 4)]   []         [(-1, 2)]\n [(-1, 4)]  []         []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]   []\n [(1, 3)]   []         [(-1, 1)]  []\n [(-1, 2)]  [(1, 1)]   []         []"
            @test repr("text/plain", PD.liealgebra_so_extpowers_standard_module_struct_const(5, 1)) ==
                  "10×5 Matrix{Vector{Tuple{Int64, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []         []         []\n [(-1, 3)]  []         [(1, 1)]   []         []\n [(-1, 4)]  []         []         [(1, 1)]   []\n [(-1, 5)]  []         []         []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]   []         []\n []         [(-1, 4)]  []         [(1, 2)]   []\n []         [(-1, 5)]  []         []         [(1, 2)]\n []         []         [(-1, 4)]  [(1, 3)]   []\n []         []         [(-1, 5)]  []         [(1, 3)]\n []         []         []         [(-1, 5)]  [(1, 4)]"
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
