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

@testset ExtendedTestSet "All LieAlgebraStructConsts.jl tests" begin
    @testset "liealgebra_so" begin
        @testset "lie algebra; n = $n" for n in 1:6
            R = QQ
            if n % 2 == 1
                l = div(n - 1, 2)
                dimL = 2 * l^2 + l
            else
                l = div(n, 2)
                dimL = 2 * l^2 - l
            end

            @test dimL == length(liealgebra_so_basis(n, R))
            @test dimL == length(liealgebra_so_symbols(n))
            @test (dimL, dimL) == size(liealgebra_so_struct_const(n, R))
        end

        @testset "symmpowers std module; n = $n, e = $e" for n in 1:5, e in 1:3
            R = QQ
            dimL = length(liealgebra_so_basis(n, R))
            dimV = binomial(n + e - 1, e)
            @test dimV == length(liealgebra_so_symmpowers_standard_module_symbols(n, e))
            @test (dimL, dimV) == size(liealgebra_so_symmpowers_standard_module_struct_const(n, e, R))
        end

        @testset "extpowers std module; n = $n, e = $e" for n in 1:5, e in 1:min(4, n)
            R = QQ
            dimL = length(liealgebra_so_basis(n, R))
            dimV = binomial(n, e)
            @test dimV == length(liealgebra_so_extpowers_standard_module_symbols(n, e))
            @test (dimL, dimV) == size(liealgebra_so_extpowers_standard_module_struct_const(n, e, R))
        end

        @testset "correctness regression" begin
            R = QQ
            @test repr("text/plain", PD.liealgebra_so_struct_const(3, R)) ==
                  "3×3 Matrix{Vector{Tuple{fmpq, Int64}}}:\n []         [(-1, 3)]  [(1, 2)]\n [(1, 3)]   []         [(-1, 1)]\n [(-1, 2)]  [(1, 1)]   []"
            @test repr("text/plain", PD.liealgebra_so_struct_const(4, R)) ==
                  "6×6 Matrix{Vector{Tuple{fmpq, Int64}}}:\n []         [(-1, 4)]  [(-1, 5)]  [(1, 2)]   [(1, 3)]   []\n [(1, 4)]   []         [(-1, 6)]  [(-1, 1)]  []         [(1, 3)]\n [(1, 5)]   [(1, 6)]   []         []         [(-1, 1)]  [(-1, 2)]\n [(-1, 2)]  [(1, 1)]   []         []         [(-1, 6)]  [(1, 5)]\n [(-1, 3)]  []         [(1, 1)]   [(1, 6)]   []         [(-1, 4)]\n []         [(-1, 3)]  [(1, 2)]   [(-1, 5)]  [(1, 4)]   []"
            @test repr("text/plain", PD.liealgebra_so_symmpowers_standard_module_struct_const(3, 1, R)) ==
                  "3×3 Matrix{Vector{Tuple{fmpq, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []\n [(-1, 3)]  []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]"
            @test repr("text/plain", PD.liealgebra_so_symmpowers_standard_module_struct_const(3, 2, R)) ==
                  "3×6 Matrix{Vector{Tuple{fmpq, Int64}}}:\n [(-2, 2)]  [(1, 1), (-1, 4)]  [(-1, 5)]          [(2, 2)]   [(1, 3)]           []\n [(-2, 3)]  [(-1, 5)]          [(1, 1), (-1, 6)]  []         [(1, 2)]           [(2, 3)]\n []         [(-1, 3)]          [(1, 2)]           [(-2, 5)]  [(1, 4), (-1, 6)]  [(2, 5)]"
            @test repr("text/plain", PD.liealgebra_so_symmpowers_standard_module_struct_const(3, 3, R)) ==
                  "3×10 Matrix{Vector{Tuple{fmpq, Int64}}}:\n [(-3, 2)]  [(1, 1), (-2, 4)]  [(-2, 5)]          [(2, 2), (-1, 7)]  [(1, 3), (-1, 8)]  [(-1, 9)]           [(3, 4)]   [(2, 5)]           [(1, 6)]            []\n [(-3, 3)]  [(-2, 5)]          [(1, 1), (-2, 6)]  [(-1, 8)]          [(1, 2), (-1, 9)]  [(2, 3), (-1, 10)]  []         [(1, 4)]           [(2, 5)]            [(3, 6)]\n []         [(-1, 3)]          [(1, 2)]           [(-2, 5)]          [(1, 4), (-1, 6)]  [(2, 5)]            [(-3, 8)]  [(1, 7), (-2, 9)]  [(2, 8), (-1, 10)]  [(3, 9)]"
            @test repr("text/plain", PD.liealgebra_so_symmpowers_standard_module_struct_const(4, 1, R)) ==
                  "6×4 Matrix{Vector{Tuple{fmpq, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []         []\n [(-1, 3)]  []         [(1, 1)]   []\n [(-1, 4)]  []         []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]   []\n []         [(-1, 4)]  []         [(1, 2)]\n []         []         [(-1, 4)]  [(1, 3)]"
            @test repr("text/plain", PD.liealgebra_so_extpowers_standard_module_struct_const(3, 1, R)) ==
                  "3×3 Matrix{Vector{Tuple{fmpq, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []\n [(-1, 3)]  []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]"
            @test repr("text/plain", PD.liealgebra_so_extpowers_standard_module_struct_const(3, 2, R)) ==
                  "3×3 Matrix{Vector{Tuple{fmpq, Int64}}}:\n []         [(-1, 3)]  [(1, 2)]\n [(1, 3)]   []         [(-1, 1)]\n [(-1, 2)]  [(1, 1)]   []"
            @test repr("text/plain", PD.liealgebra_so_extpowers_standard_module_struct_const(4, 1, R)) ==
                  "6×4 Matrix{Vector{Tuple{fmpq, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []         []\n [(-1, 3)]  []         [(1, 1)]   []\n [(-1, 4)]  []         []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]   []\n []         [(-1, 4)]  []         [(1, 2)]\n []         []         [(-1, 4)]  [(1, 3)]"
            @test repr("text/plain", PD.liealgebra_so_extpowers_standard_module_struct_const(4, 2, R)) ==
                  "6×6 Matrix{Vector{Tuple{fmpq, Int64}}}:\n []         [(-1, 4)]  [(-1, 5)]  [(1, 2)]   [(1, 3)]   []\n [(1, 4)]   []         [(-1, 6)]  [(-1, 1)]  []         [(1, 3)]\n [(1, 5)]   [(1, 6)]   []         []         [(-1, 1)]  [(-1, 2)]\n [(-1, 2)]  [(1, 1)]   []         []         [(-1, 6)]  [(1, 5)]\n [(-1, 3)]  []         [(1, 1)]   [(1, 6)]   []         [(-1, 4)]\n []         [(-1, 3)]  [(1, 2)]   [(-1, 5)]  [(1, 4)]   []"
            @test repr("text/plain", PD.liealgebra_so_extpowers_standard_module_struct_const(4, 3, R)) ==
                  "6×4 Matrix{Vector{Tuple{fmpq, Int64}}}:\n []         []         [(-1, 4)]  [(1, 3)]\n []         [(1, 4)]   []         [(-1, 2)]\n [(-1, 4)]  []         []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]   []\n [(1, 3)]   []         [(-1, 1)]  []\n [(-1, 2)]  [(1, 1)]   []         []"
            @test repr("text/plain", PD.liealgebra_so_extpowers_standard_module_struct_const(5, 1, R)) ==
                  "10×5 Matrix{Vector{Tuple{fmpq, Int64}}}:\n [(-1, 2)]  [(1, 1)]   []         []         []\n [(-1, 3)]  []         [(1, 1)]   []         []\n [(-1, 4)]  []         []         [(1, 1)]   []\n [(-1, 5)]  []         []         []         [(1, 1)]\n []         [(-1, 3)]  [(1, 2)]   []         []\n []         [(-1, 4)]  []         [(1, 2)]   []\n []         [(-1, 5)]  []         []         [(1, 2)]\n []         []         [(-1, 4)]  [(1, 3)]   []\n []         []         [(-1, 5)]  []         [(1, 3)]\n []         []         []         [(-1, 5)]  [(1, 4)]"
        end

    end


end
