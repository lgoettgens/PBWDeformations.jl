@testset "ActingGroup.jl tests" begin

    function is_trivial_hom(h, G)
        @assert G === domain(h)
        return all(isone ∘ h, G)
    end

    function is_usual_sign(h, G)
        @assert G === domain(h)
        return all(g -> isone(h(g)) == isone(sign(g)), G)
    end

    @testset "acting_group_with_sgn; $L" for L in (
        special_orthogonal_lie_algebra(QQ, 3, identity_matrix(QQ, 3)),
        general_linear_lie_algebra(QQ, 2),
    )
        Vstd = standard_module(L)

        @testset "standard module" begin
            V = Vstd
            G, h = acting_group_with_sgn(V)
            @test G == symmetric_group(1)
            @test is_trivial_hom(h, G)
        end

        @testset "dual standard module" begin
            V = dual(Vstd)
            G, h = acting_group_with_sgn(V)
            @test G == symmetric_group(1)
            @test is_trivial_hom(h, G)
        end

        @testset "exterior power" begin
            @testset "exterior power ($k-th)" for k in 1:4
                V = exterior_power_obj(Vstd, k)
                G, h = acting_group_with_sgn(V)
                @test G == symmetric_group(k)
                @test is_usual_sign(h, G)
            end
        end

        @testset "symmetric power" begin
            @testset "symmetric power ($k-th)" for k in 1:4
                V = symmetric_power_obj(Vstd, k)
                G, h = acting_group_with_sgn(V)
                @test G == symmetric_group(k)
                @test is_trivial_hom(h, G)
            end
        end

        @testset "tensor power" begin
            @testset "tensor power ($k-th)" for k in 1:4
                V = tensor_power_obj(Vstd, k)
                G, h = acting_group_with_sgn(V)
                @test G == permutation_group(k, PermGroupElem[])
                @test is_trivial_hom(h, G)
            end
        end

        @testset "tensor product" begin
            V = tensor_product(exterior_power_obj(Vstd, 3), symmetric_power_obj(Vstd, 3))
            G, h = acting_group_with_sgn(V)
            @test G == permutation_group(6,[cperm(G, 1:2),cperm(G, 1:3),cperm(G, 3 .+ (1:2)),cperm(G, 3 .+ (1:3))])
            @test !isone(h(cperm(G, 1:2)))
            @test isone(h(cperm(G, 3 .+ (1:2))))
        end

        @testset "exterior power of exterior power" begin
            @testset "exterior power ($l-th) of exterior power ($k-th)" for k in 1:5, l in 1:5
                V = exterior_power_obj(exterior_power_obj(Vstd, k), l)
                G, h = acting_group_with_sgn(V)
                @test G == permutation_group(
                    k * l,
                    [
                        cperm(1:k),
                        cperm(1:min(k, 2)),
                        cperm([[(i - 1) * k + j for i in 1:l] for j in 1:k]),
                        cperm([[(i - 1) * k + j for i in 1:min(l, 2)] for j in 1:k]),
                    ],
                )
                if k > 1
                    for i in 1:l
                        @test !isone(h(cperm(G, (i-1)*k .+ (1:2))))
                        @test isone(h(cperm(G, (i-1)*k .+ (1:k)))) == isodd(k)
                    end
                else
                    @test is_usual_sign(h, G)
                end
                if l > 1
                    @test !isone(h(cperm([[(i - 1) * k + j for i in 1:2] for j in 1:k])))
                    @test isone(h(cperm([[(i - 1) * k + j for i in 1:l] for j in 1:k]))) == isodd(l)
                else
                    @test is_usual_sign(h, G)
                end
            end
        end

        @testset "symmetric power of symmetric power" begin
            @testset "symmetric power ($l-th) of symmetric power ($k-th)" for (k,l) in [(k,l) for k in 1:4 for l in 1:4 if k*(l+1) < 12]
                V = symmetric_power_obj(symmetric_power_obj(Vstd, k), l)
                G, h = acting_group_with_sgn(V)
                @test G == permutation_group(
                    k * l,
                    [
                        cperm(1:k),
                        cperm(1:min(k, 2)),
                        cperm([[(i - 1) * k + j for i in 1:l] for j in 1:k]),
                        cperm([[(i - 1) * k + j for i in 1:min(l, 2)] for j in 1:k]),
                    ],
                )
                @test is_trivial_hom(h, G)
            end
        end

        @testset "tensor power of tensor power" begin
            @testset "tensor power ($l-th) of tensor power ($k-th)" for (k,l) in [(k,l) for k in 1:4 for l in 1:4 if k*l < 8]
                V = tensor_power_obj(tensor_power_obj(Vstd, k), l)
                G, h = acting_group_with_sgn(V)
                @test G == permutation_group(k*l, PermGroupElem[])
                @test is_trivial_hom(h, G)
            end
        end

        @testset "symmetric power of exterior power" begin
            @testset "symmetric power ($l-th) of exterior power ($k-th)" for k in 1:5, l in 1:5
                if k > dim(Vstd) # TODO: remove once https://github.com/oscar-system/Oscar.jl/pull/4878 is available
                    continue
                end
                V = symmetric_power_obj(exterior_power_obj(Vstd, k), l)
                G, h = acting_group_with_sgn(V)
                @test G == permutation_group(
                    k * l,
                    [
                        cperm(1:k),
                        cperm(1:min(k, 2)),
                        cperm([[(i - 1) * k + j for i in 1:l] for j in 1:k]),
                        cperm([[(i - 1) * k + j for i in 1:min(l, 2)] for j in 1:k]),
                    ],
                )
                if k > 1
                    for i in 1:l
                        @test !isone(h(cperm(G, (i-1)*k .+ (1:2))))
                        @test isone(h(cperm(G, (i-1)*k .+ (1:k)))) == isodd(k)
                    end
                else
                    @test is_trivial_hom(h, G)
                end
                if l > 1
                    @test isone(h(cperm([[(i - 1) * k + j for i in 1:2] for j in 1:k])))
                    @test isone(h(cperm([[(i - 1) * k + j for i in 1:l] for j in 1:k])))
                else
                    @test is_usual_sign(h, G)
                end
            end
        end

        @testset "exterior power of symmetric power" begin
            @testset "exterior power ($l-th) of symmetric power ($k-th)" for (k,l) in [(k,l) for k in 1:5 for l in 1:4 if k*l < 12]
                V = exterior_power_obj(symmetric_power_obj(Vstd, k), l)
                G, h = acting_group_with_sgn(V)
                @test G == permutation_group(
                    k * l,
                    [
                        cperm(1:k),
                        cperm(1:min(k, 2)),
                        cperm([[(i - 1) * k + j for i in 1:l] for j in 1:k]),
                        cperm([[(i - 1) * k + j for i in 1:min(l, 2)] for j in 1:k]),
                    ],
                )
                if k > 1
                    for i in 1:l
                        @test isone(h(cperm(G, (i-1)*k .+ (1:2))))
                        @test isone(h(cperm(G, (i-1)*k .+ (1:k))))
                    end
                else
                    @test is_usual_sign(h, G)
                end
                if l > 1
                    @test !isone(h(cperm([[(i - 1) * k + j for i in 1:2] for j in 1:k])))
                    @test isone(h(cperm([[(i - 1) * k + j for i in 1:l] for j in 1:k]))) == isodd(l)
                else
                    @test is_trivial_hom(h, G)
                end
            end
        end

        @testset "tensor power of exterior power" begin
            @testset "tensor power ($l-th) of exterior power ($k-th)" for k in 1:3, l in 1:5
                V = tensor_power_obj(exterior_power_obj(Vstd, k), l)
                G, h = acting_group_with_sgn(V)
                @test G == permutation_group(
                    k * l,
                    [
                        (cperm((i - 1) * k .+ (1:k)) for i in 1:l)...,
                        (cperm((i - 1) * k .+ (1:min(k, 2))) for i in 1:l)...,
                    ],
                )
                if k > 1
                    for i in 1:l
                        @test !isone(h(cperm(G, (i-1)*k .+ (1:2))))
                        @test isone(h(cperm(G, (i-1)*k .+ (1:k)))) == isodd(k)
                    end
                else
                    @test is_trivial_hom(h, G)
                end
            end
        end

        @testset "exterior power of tensor power" begin
            @testset "exterior power ($l-th) of tensor power ($k-th)" for (k,l) in [(k,l) for k in 1:3 for l in 1:4 if (k+1)*l < 12]
                V = exterior_power_obj(tensor_power_obj(Vstd, k), l)
                G, h = acting_group_with_sgn(V)
                @test G == permutation_group(
                    k * l,
                    [
                        cperm([[(i - 1) * k + j for i in 1:l] for j in 1:k]),
                        cperm([[(i - 1) * k + j for i in 1:min(l, 2)] for j in 1:k]),
                    ],
                )
                if l > 1
                    @test !isone(h(cperm([[(i - 1) * k + j for i in 1:2] for j in 1:k])))
                    @test isone(h(cperm([[(i - 1) * k + j for i in 1:l] for j in 1:k]))) == isodd(l)
                else
                    @test is_trivial_hom(h, G)
                end
            end
        end

        @testset "tensor power of symmetric power" begin
            @testset "tensor power ($l-th) of symmetric power ($k-th)" for (k,l) in [(k,l) for k in 1:3 for l in 1:4 if k*(l+1) < 12]
                V = tensor_power_obj(symmetric_power_obj(Vstd, k), l)
                G, h = acting_group_with_sgn(V)
                @test G == permutation_group(
                    k * l,
                    [
                        (cperm((i - 1) * k .+ (1:k)) for i in 1:l)...,
                        (cperm((i - 1) * k .+ (1:min(k, 2))) for i in 1:l)...,
                    ],
                )
                if k > 1
                    for i in 1:l
                        @test isone(h(cperm(G, (i-1)*k .+ (1:2))))
                        @test isone(h(cperm(G, (i-1)*k .+ (1:k))))
                    end
                else
                    @test is_trivial_hom(h, G)
                end
            end
        end

        @testset "symmetric power of tensor power" begin
            @testset "symmetric power ($l-th) of tensor power ($k-th)" for (k,l) in [(k,l) for k in 1:3 for l in 1:4 if (k+1)*l < 12]
                V = symmetric_power_obj(tensor_power_obj(Vstd, k), l)
                G, h = acting_group_with_sgn(V)
                @test G == permutation_group(
                    k * l,
                    [
                        cperm([[(i - 1) * k + j for i in 1:l] for j in 1:k]),
                        cperm([[(i - 1) * k + j for i in 1:min(l, 2)] for j in 1:k]),
                    ],
                )
                if l > 1
                    @test isone(h(cperm([[(i - 1) * k + j for i in 1:2] for j in 1:k])))
                    @test isone(h(cperm([[(i - 1) * k + j for i in 1:l] for j in 1:k])))
                else
                    @test is_trivial_hom(h, G)
                end
            end
        end
    end
end
