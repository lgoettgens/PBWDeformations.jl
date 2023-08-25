@testset "ModuleSimpleStructure.jl tests" begin
    @testset for L in [general_linear_lie_algebra(QQ, 3), special_orthogonal_lie_algebra(QQ, 3)]
        stdV = standard_module(L)

        @testset "Duality, standard" begin
            V = dual(stdV)
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == V
            @test sprint(show, basis(V)) == "LieAlgebraModuleElem{QQFieldElem}[v_1*, v_2*, v_3*]"
            @test sprint(show, h.(basis(V))) == "LieAlgebraModuleElem{QQFieldElem}[v_1*, v_2*, v_3*]"
        end

        @testset "Duality, duality" begin
            V = dual(dual(stdV))
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == stdV
            @test sprint(show, basis(V)) == "LieAlgebraModuleElem{QQFieldElem}[(v_1*)*, (v_2*)*, (v_3*)*]"
            @test sprint(show, h.(basis(V))) == "LieAlgebraModuleElem{QQFieldElem}[v_1, v_2, v_3]"
        end

        @testset "Duality, direct sum" begin
            V = dual(direct_sum(stdV, stdV, dual(stdV)))
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == direct_sum(dual(stdV), dual(stdV), stdV)
            @test sprint(show, basis(V)) ==
                  "LieAlgebraModuleElem{QQFieldElem}[(v_1^(1))*, (v_2^(1))*, (v_3^(1))*, (v_1^(2))*, (v_2^(2))*, (v_3^(2))*, ((v_1*)^(3))*, ((v_2*)^(3))*, ((v_3*)^(3))*]"
            @test sprint(show, h.(basis(V))) ==
                  "LieAlgebraModuleElem{QQFieldElem}[(v_1*)^(1), (v_2*)^(1), (v_3*)^(1), (v_1*)^(2), (v_2*)^(2), (v_3*)^(2), v_1^(3), v_2^(3), v_3^(3)]"
        end

        @testset "Duality, tensor product" begin
            V = dual(tensor_product(stdV, dual(stdV)))
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == tensor_product(dual(stdV), stdV)
            @test sprint(show, basis(V)) ==
                  "LieAlgebraModuleElem{QQFieldElem}[(v_1 ⊗ (v_1*))*, (v_1 ⊗ (v_2*))*, (v_1 ⊗ (v_3*))*, (v_2 ⊗ (v_1*))*, (v_2 ⊗ (v_2*))*, (v_2 ⊗ (v_3*))*, (v_3 ⊗ (v_1*))*, (v_3 ⊗ (v_2*))*, (v_3 ⊗ (v_3*))*]"
            @test sprint(show, h.(basis(V))) ==
                  "LieAlgebraModuleElem{QQFieldElem}[(v_1*) ⊗ v_1, (v_1*) ⊗ v_2, (v_1*) ⊗ v_3, (v_2*) ⊗ v_1, (v_2*) ⊗ v_2, (v_2*) ⊗ v_3, (v_3*) ⊗ v_1, (v_3*) ⊗ v_2, (v_3*) ⊗ v_3]"
        end

        @testset "Duality, exterior power, k = $k" for k in 2:3
            V = dual(exterior_power(stdV, k))
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == exterior_power(dual(stdV), k)
            if k == 2
                @test sprint(show, basis(V)) ==
                      "LieAlgebraModuleElem{QQFieldElem}[(v_1 ∧ v_2)*, (v_1 ∧ v_3)*, (v_2 ∧ v_3)*]"
                @test sprint(show, h.(basis(V))) ==
                      "LieAlgebraModuleElem{QQFieldElem}[(v_1*) ∧ (v_2*), (v_1*) ∧ (v_3*), (v_2*) ∧ (v_3*)]"
            end
        end

        @testset "Duality, symmetric power, k = $k" for k in 2:5
            V = dual(symmetric_power(stdV, k))
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == symmetric_power(dual(stdV), k)
            if k == 2
                @test sprint(show, basis(V)) ==
                      "LieAlgebraModuleElem{QQFieldElem}[(v_1^2)*, (v_1*v_2)*, (v_1*v_3)*, (v_2^2)*, (v_2*v_3)*, (v_3^2)*]"
                @test sprint(show, h.(basis(V))) ==
                      "LieAlgebraModuleElem{QQFieldElem}[(v_1*)^2, 2*(v_1*)*(v_2*), 2*(v_1*)*(v_3*), (v_2*)^2, 2*(v_2*)*(v_3*), (v_3*)^2]"
            end
        end

        @testset "Duality, tensor power, k = $k" for k in 2:5
            V = dual(tensor_power(stdV, k))
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == tensor_power(dual(stdV), k)
            if k == 2
                @test sprint(show, basis(V)) ==
                      "LieAlgebraModuleElem{QQFieldElem}[(v_1 ⊗ v_1)*, (v_1 ⊗ v_2)*, (v_1 ⊗ v_3)*, (v_2 ⊗ v_1)*, (v_2 ⊗ v_2)*, (v_2 ⊗ v_3)*, (v_3 ⊗ v_1)*, (v_3 ⊗ v_2)*, (v_3 ⊗ v_3)*]"
                @test sprint(show, h.(basis(V))) ==
                      "LieAlgebraModuleElem{QQFieldElem}[(v_1*) ⊗ (v_1*), (v_1*) ⊗ (v_2*), (v_1*) ⊗ (v_3*), (v_2*) ⊗ (v_1*), (v_2*) ⊗ (v_2*), (v_2*) ⊗ (v_3*), (v_3*) ⊗ (v_1*), (v_3*) ⊗ (v_2*), (v_3*) ⊗ (v_3*)]"
            end
        end

        @testset "Direct sum: one summand" begin
            V = direct_sum(exterior_power(stdV, 2))
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == exterior_power(stdV, 2)
            @test sprint(show, basis(V)) == "LieAlgebraModuleElem{QQFieldElem}[v_1 ∧ v_2, v_1 ∧ v_3, v_2 ∧ v_3]"
            @test sprint(show, h.(basis(V))) == "LieAlgebraModuleElem{QQFieldElem}[v_1 ∧ v_2, v_1 ∧ v_3, v_2 ∧ v_3]"
        end

        @testset "Direct sum: nested" begin
            V = direct_sum(direct_sum(stdV, exterior_power(stdV, 2)), dual(stdV))
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == direct_sum(stdV, exterior_power(stdV, 2), dual(stdV))
            @test sprint(show, basis(V)) ==
                  "LieAlgebraModuleElem{QQFieldElem}[(v_1^(1))^(1), (v_2^(1))^(1), (v_3^(1))^(1), ((v_1 ∧ v_2)^(2))^(1), ((v_1 ∧ v_3)^(2))^(1), ((v_2 ∧ v_3)^(2))^(1), (v_1*)^(2), (v_2*)^(2), (v_3*)^(2)]"
            @test sprint(show, h.(basis(V))) ==
                  "LieAlgebraModuleElem{QQFieldElem}[v_1^(1), v_2^(1), v_3^(1), (v_1 ∧ v_2)^(2), (v_1 ∧ v_3)^(2), (v_2 ∧ v_3)^(2), (v_1*)^(3), (v_2*)^(3), (v_3*)^(3)]"
        end

        @testset "Tensor product: one factor" begin
            V = tensor_product(exterior_power(stdV, 2))
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == exterior_power(stdV, 2)
            @test sprint(show, basis(V)) == "LieAlgebraModuleElem{QQFieldElem}[v_1 ∧ v_2, v_1 ∧ v_3, v_2 ∧ v_3]"
            @test sprint(show, h.(basis(V))) == "LieAlgebraModuleElem{QQFieldElem}[v_1 ∧ v_2, v_1 ∧ v_3, v_2 ∧ v_3]"
        end

        @testset "Tensor product: nested" begin
            V = tensor_product(tensor_product(stdV, exterior_power(stdV, 3)), dual(stdV))
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == tensor_product(stdV, exterior_power(stdV, 3), dual(stdV))
            @test sprint(show, basis(V)) ==
                  "LieAlgebraModuleElem{QQFieldElem}[(v_1 ⊗ (v_1 ∧ v_2 ∧ v_3)) ⊗ (v_1*), (v_1 ⊗ (v_1 ∧ v_2 ∧ v_3)) ⊗ (v_2*), (v_1 ⊗ (v_1 ∧ v_2 ∧ v_3)) ⊗ (v_3*), (v_2 ⊗ (v_1 ∧ v_2 ∧ v_3)) ⊗ (v_1*), (v_2 ⊗ (v_1 ∧ v_2 ∧ v_3)) ⊗ (v_2*), (v_2 ⊗ (v_1 ∧ v_2 ∧ v_3)) ⊗ (v_3*), (v_3 ⊗ (v_1 ∧ v_2 ∧ v_3)) ⊗ (v_1*), (v_3 ⊗ (v_1 ∧ v_2 ∧ v_3)) ⊗ (v_2*), (v_3 ⊗ (v_1 ∧ v_2 ∧ v_3)) ⊗ (v_3*)]"
            @test sprint(show, h.(basis(V))) ==
                  "LieAlgebraModuleElem{QQFieldElem}[v_1 ⊗ (v_1 ∧ v_2 ∧ v_3) ⊗ (v_1*), v_1 ⊗ (v_1 ∧ v_2 ∧ v_3) ⊗ (v_2*), v_1 ⊗ (v_1 ∧ v_2 ∧ v_3) ⊗ (v_3*), v_2 ⊗ (v_1 ∧ v_2 ∧ v_3) ⊗ (v_1*), v_2 ⊗ (v_1 ∧ v_2 ∧ v_3) ⊗ (v_2*), v_2 ⊗ (v_1 ∧ v_2 ∧ v_3) ⊗ (v_3*), v_3 ⊗ (v_1 ∧ v_2 ∧ v_3) ⊗ (v_1*), v_3 ⊗ (v_1 ∧ v_2 ∧ v_3) ⊗ (v_2*), v_3 ⊗ (v_1 ∧ v_2 ∧ v_3) ⊗ (v_3*)]"
        end

        @testset "Tensor product, direct sum" begin
            V = tensor_product(direct_sum(stdV, exterior_power(stdV, 3)), direct_sum(dual(stdV), stdV))
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == direct_sum(
                tensor_product(stdV, dual(stdV)),
                tensor_product(stdV, stdV),
                tensor_product(exterior_power(stdV, 3), dual(stdV)),
                tensor_product(exterior_power(stdV, 3), stdV),
            )
            @test sprint(show, basis(V)) ==
                  "LieAlgebraModuleElem{QQFieldElem}[(v_1^(1)) ⊗ ((v_1*)^(1)), (v_1^(1)) ⊗ ((v_2*)^(1)), (v_1^(1)) ⊗ ((v_3*)^(1)), (v_1^(1)) ⊗ (v_1^(2)), (v_1^(1)) ⊗ (v_2^(2)), (v_1^(1)) ⊗ (v_3^(2)), (v_2^(1)) ⊗ ((v_1*)^(1)), (v_2^(1)) ⊗ ((v_2*)^(1)), (v_2^(1)) ⊗ ((v_3*)^(1)), (v_2^(1)) ⊗ (v_1^(2)), (v_2^(1)) ⊗ (v_2^(2)), (v_2^(1)) ⊗ (v_3^(2)), (v_3^(1)) ⊗ ((v_1*)^(1)), (v_3^(1)) ⊗ ((v_2*)^(1)), (v_3^(1)) ⊗ ((v_3*)^(1)), (v_3^(1)) ⊗ (v_1^(2)), (v_3^(1)) ⊗ (v_2^(2)), (v_3^(1)) ⊗ (v_3^(2)), ((v_1 ∧ v_2 ∧ v_3)^(2)) ⊗ ((v_1*)^(1)), ((v_1 ∧ v_2 ∧ v_3)^(2)) ⊗ ((v_2*)^(1)), ((v_1 ∧ v_2 ∧ v_3)^(2)) ⊗ ((v_3*)^(1)), ((v_1 ∧ v_2 ∧ v_3)^(2)) ⊗ (v_1^(2)), ((v_1 ∧ v_2 ∧ v_3)^(2)) ⊗ (v_2^(2)), ((v_1 ∧ v_2 ∧ v_3)^(2)) ⊗ (v_3^(2))]"
            @test sprint(show, h.(basis(V))) ==
                  "LieAlgebraModuleElem{QQFieldElem}[(v_1 ⊗ (v_1*))^(1), (v_1 ⊗ (v_2*))^(1), (v_1 ⊗ (v_3*))^(1), (v_1 ⊗ v_1)^(2), (v_1 ⊗ v_2)^(2), (v_1 ⊗ v_3)^(2), (v_2 ⊗ (v_1*))^(1), (v_2 ⊗ (v_2*))^(1), (v_2 ⊗ (v_3*))^(1), (v_2 ⊗ v_1)^(2), (v_2 ⊗ v_2)^(2), (v_2 ⊗ v_3)^(2), (v_3 ⊗ (v_1*))^(1), (v_3 ⊗ (v_2*))^(1), (v_3 ⊗ (v_3*))^(1), (v_3 ⊗ v_1)^(2), (v_3 ⊗ v_2)^(2), (v_3 ⊗ v_3)^(2), ((v_1 ∧ v_2 ∧ v_3) ⊗ (v_1*))^(3), ((v_1 ∧ v_2 ∧ v_3) ⊗ (v_2*))^(3), ((v_1 ∧ v_2 ∧ v_3) ⊗ (v_3*))^(3), ((v_1 ∧ v_2 ∧ v_3) ⊗ v_1)^(4), ((v_1 ∧ v_2 ∧ v_3) ⊗ v_2)^(4), ((v_1 ∧ v_2 ∧ v_3) ⊗ v_3)^(4)]"


            V1 = exterior_power(stdV, 2)
            V2 = dual(stdV)
            V3 = symmetric_power(stdV, 2)
            V4 = dual(stdV)
            V5 = stdV
            V6 = exterior_power(stdV, 3)
            V = tensor_product(V1, direct_sum(V2, V3, V4), direct_sum(V5, V6))
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == direct_sum(
                tensor_product(V1, V2, V5),
                tensor_product(V1, V2, V6),
                tensor_product(V1, V3, V5),
                tensor_product(V1, V3, V6),
                tensor_product(V1, V4, V5),
                tensor_product(V1, V4, V6),
            )
        end

        @testset "Exterior power: 1st power" begin
            V1 = stdV
            V = exterior_power(V1, 1)
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == V1
        end

        @testset "Exterior power, direct sum" begin
            V1 = stdV
            V2 = dual(stdV)
            V3 = exterior_power(stdV, 2)
            V = exterior_power(direct_sum(V1, V2, V3), 3)
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == direct_sum(
                exterior_power(V1, 3),
                tensor_product(exterior_power(V1, 2), V2),
                tensor_product(exterior_power(V1, 2), V3),
                tensor_product(V1, exterior_power(V2, 2)),
                tensor_product(V1, V2, V3),
                tensor_product(V1, exterior_power(V3, 2)),
                exterior_power(V2, 3),
                tensor_product(exterior_power(V2, 2), V3),
                tensor_product(V2, exterior_power(V3, 2)),
                exterior_power(V3, 3),
            )
        end

        @testset "Symmetric power: 1st power" begin
            V1 = stdV
            V = symmetric_power(V1, 1)
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == V1
        end

        @testset "Symmetric power, direct sum" begin
            V1 = stdV
            V2 = dual(stdV)
            V3 = exterior_power(stdV, 2)
            V = symmetric_power(direct_sum(V1, V2, V3), 3)
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == direct_sum(
                symmetric_power(V1, 3),
                tensor_product(symmetric_power(V1, 2), V2),
                tensor_product(symmetric_power(V1, 2), V3),
                tensor_product(V1, symmetric_power(V2, 2)),
                tensor_product(V1, V2, V3),
                tensor_product(V1, symmetric_power(V3, 2)),
                symmetric_power(V2, 3),
                tensor_product(symmetric_power(V2, 2), V3),
                tensor_product(V2, symmetric_power(V3, 2)),
                symmetric_power(V3, 3),
            )
        end

        @testset "Tensor power: 1st power" begin
            V1 = stdV
            V = tensor_power(V1, 1)
            W, h = isomorphic_module_with_simple_structure(V)
            @test is_isomorphism(h)
            @test W == V1
        end

        @testset "Tensor power, direct sum" begin
            V1 = stdV
            V2 = exterior_power(dual(stdV), 3)
            V3 = exterior_power(stdV, 2)
            V = tensor_power(direct_sum(V1, V2, V3), 3)
            W, h = isomorphic_module_with_simple_structure(V; check=false) # otherwise needs 3 minutes
            @test is_isomorphism(h)
            @test W == direct_sum(
                tensor_power(V1, 3),                        # 111
                tensor_product(tensor_power(V1, 2), V2),    # 112
                tensor_product(tensor_power(V1, 2), V3),    # 113
                tensor_product(tensor_power(V1, 2), V2),    # 121
                tensor_product(V1, tensor_power(V2, 2)),    # 122
                tensor_product(V1, V2, V3),                 # 123
                tensor_product(tensor_power(V1, 2), V3),    # 131
                tensor_product(V1, V2, V3),                 # 132
                tensor_product(V1, tensor_power(V3, 2)),    # 133
                tensor_product(tensor_power(V1, 2), V2),    # 211
                tensor_product(V1, tensor_power(V2, 2)),    # 212
                tensor_product(V1, V2, V3),                 # 213
                tensor_product(V1, tensor_power(V2, 2)),    # 221
                tensor_power(V2, 3),                        # 222
                tensor_product(tensor_power(V2, 2), V3),    # 223
                tensor_product(V1, V2, V3),                 # 231
                tensor_product(tensor_power(V2, 2), V3),    # 232
                tensor_product(V2, tensor_power(V3, 2)),    # 233
                tensor_product(tensor_power(V1, 2), V3),    # 311
                tensor_product(V1, V2, V3),                 # 312
                tensor_product(V1, tensor_power(V3, 2)),    # 313
                tensor_product(V1, V2, V3),                 # 321
                tensor_product(tensor_power(V2, 2), V3),    # 322
                tensor_product(V2, tensor_power(V3, 2)),    # 323
                tensor_product(V1, tensor_power(V3, 2)),    # 331
                tensor_product(V2, tensor_power(V3, 2)),    # 332
                tensor_power(V3, 3),                        # 333
            )
        end

        @testset "Complex cases 1" begin
            V1 = stdV
            V2 = dual(V1)
            V3 = direct_sum(V1, V2)

            V_1 = exterior_power(V3, 1)
            W_1, h = isomorphic_module_with_simple_structure(V_1)
            @test is_isomorphism(h)
            @test W_1 == V3

            V_2 = exterior_power(V3, 2)
            W_2, h = isomorphic_module_with_simple_structure(V_2)
            @test is_isomorphism(h)
            @test W_2 == direct_sum(exterior_power(V1, 2), tensor_product(V1, V2), exterior_power(V2, 2))

            V_3 = exterior_power(V3, 3)
            W_3, h = isomorphic_module_with_simple_structure(V_3)
            @test is_isomorphism(h)
            @test W_3 == direct_sum(
                exterior_power(V1, 3),
                tensor_product(exterior_power(V1, 2), V2),
                tensor_product(V1, exterior_power(V2, 2)),
                exterior_power(V2, 3),
            )

            # TODO: handle zero-dim modules (e.g. 4th exterior power of 3-dim module)
            # V_4 = exterior_power(V3, 4)
            # W_4, h = isomorphic_module_with_simple_structure(V_4)
            # @test is_isomorphism(h)
            # @test W_4 == direct_sum(
            #     exterior_power(V1, 4),
            #     tensor_product(exterior_power(V1, 3), V2),
            #     tensor_product(exterior_power(V1, 2), exterior_power(V2, 2)),
            #     tensor_product(V1, exterior_power(V2, 3)),
            #     exterior_power(V2, 4),
            # )
        end

        @testset "Complex cases 2" begin
            V1 = symmetric_power(stdV, 3)

            Ve2 = dual(exterior_power(dual(V1), 2))
            We2, h = isomorphic_module_with_simple_structure(Ve2)
            @test is_isomorphism(h)
            @test We2 == exterior_power(V1, 2)

            Ve3 = dual(exterior_power(dual(V1), 3))
            We3, h = isomorphic_module_with_simple_structure(Ve3)
            @test is_isomorphism(h)
            @test We3 == exterior_power(V1, 3)

            Vs2 = dual(symmetric_power(dual(V1), 2))
            Ws2, h = isomorphic_module_with_simple_structure(Vs2)
            @test is_isomorphism(h)
            @test Ws2 == symmetric_power(V1, 2)

            Vt2 = dual(tensor_power(dual(V1), 2))
            Wt2, h = isomorphic_module_with_simple_structure(Vt2)
            @test is_isomorphism(h)
            @test Wt2 == tensor_power(V1, 2)
        end
    end
end
