@testset verbose=true "Database.jl tests" begin
    @testset "string_for_filename" begin
        string_for_filename = PBWDeformations.Database.string_for_filename
        string_for_filename_setup = PBWDeformations.Database.string_for_filename_setup
        string_for_filename_pbwdeforms = PBWDeformations.Database.string_for_filename_pbwdeforms

        @test string_for_filename(general_linear_lie_algebra(QQ, 3)) == "gl_3"
        @test string_for_filename(special_orthogonal_lie_algebra(QQ, 4)) == "so_4"

        L = general_linear_lie_algebra(QQ, 2)

        @test string_for_filename(standard_module(L)) == "V"
        @test string_for_filename(dual(standard_module(L))) == "DV"
        @test string_for_filename(direct_sum(standard_module(L), dual(standard_module(L)))) == "_V_+_DV_"
        @test string_for_filename(tensor_product(standard_module(L), dual(standard_module(L)))) == "_V_x_DV_"
        @test string_for_filename(exterior_power_obj(standard_module(L), 2)) == "E2_V_"
        @test string_for_filename(exterior_power_obj(dual(standard_module(L)), 3)) == "E3_DV_"
        @test string_for_filename(dual(exterior_power_obj(standard_module(L), 2))) == "D_E2_V__"
        @test string_for_filename(symmetric_power_obj(standard_module(L), 2)) == "S2_V_"
        @test string_for_filename(symmetric_power_obj(dual(standard_module(L)), 2)) == "S2_DV_"
        @test string_for_filename(dual(symmetric_power_obj(standard_module(L), 3))) == "D_S3_V__"
        @test string_for_filename(tensor_power_obj(standard_module(L), 1)) == "T1_V_"
        @test string_for_filename(tensor_power_obj(dual(standard_module(L)), 2)) == "T2_DV_"
        @test string_for_filename(dual(tensor_power_obj(standard_module(L), 2))) == "D_T2_V__"
        @test string_for_filename(exterior_power_obj(direct_sum(standard_module(L), dual(standard_module(L))), 2)) == "E2__V_+_DV__"
        @test string_for_filename(direct_sum(exterior_power_obj(standard_module(L), 2), exterior_power_obj(dual(standard_module(L)), 3))) == "_E2_V__+_E3_DV__"
        @test string_for_filename(exterior_power_obj(tensor_product(standard_module(L), dual(standard_module(L))), 2)) == "E2__V_x_DV__"
        @test string_for_filename(tensor_product(exterior_power_obj(standard_module(L), 2), exterior_power_obj(dual(standard_module(L)), 3))) == "_E2_V__x_E3_DV__"

        V = direct_sum(standard_module(L), dual(standard_module(L)))
        sp = smash_product(L, V)
        @test string_for_filename(sp) == "gl_2-_V_+_DV_"
        @test string_for_filename_setup(sp) == "setup-gl_2-_V_+_DV_"

        @test string_for_filename(ArcDiagDeformBasis, sp, [1, 2]) == "ArcDiagDeformBasis-gl_2-_V_+_DV_-1_2"
        @test string_for_filename(GlnGraphDeformBasis, sp, 3:3) == "GlnGraphDeformBasis-gl_2-_V_+_DV_-3"
        @test string_for_filename_pbwdeforms(sp, [1, 2]) == "PBWDeformations-gl_2-_V_+_DV_-1_2"
        @test string_for_filename_pbwdeforms(sp, 3:3) == "PBWDeformations-gl_2-_V_+_DV_-3"
    end

    @testset verbose=true "saving and loading" begin
        compute_and_save_instance = PBWDeformations.Database.compute_and_save_instance
        load_glngraph_deform_basis = PBWDeformations.Database.load_glngraph_deform_basis
        load_glngraph_deform_bases = PBWDeformations.Database.load_glngraph_deform_bases
        load_pbwdeformations = PBWDeformations.Database.load_pbwdeformations

        mktempdir() do db
            L = general_linear_lie_algebra(QQ, 2)
            V = tensor_product(standard_module(L), dual(standard_module(L)))
            sp = smash_product(L, V)

            compute_and_save_instance(db, sp, 4)

            for reset_inbetween in [false, true]
                @testset "loading GlnGraphDeformBasis; reset_inbetween=$(reset_inbetween)" begin
                    reset_inbetween && Oscar.reset_global_serializer_state()
                    b = load_glngraph_deform_basis(db, sp, 2:2)
                    @test b isa GlnGraphDeformBasis
                    @test b.degs == 2:2
                    reset_inbetween || @test b.sp == sp
                    b_new = GlnGraphDeformBasis(sp, 2:2)
                    @test length(b) == length(b_new)
                    reset_inbetween || @test collect(b) == collect(b_new)

                    reset_inbetween && Oscar.reset_global_serializer_state()
                    bs = load_glngraph_deform_bases(db, sp, [0:3, 1:1])
                    @test all(b -> b isa GlnGraphDeformBasis, bs)
                    @test length(bs) == 2
                    @test map(b -> b.degs, bs) == [0:3, 1:1]
                    reset_inbetween || @test all(b -> b.sp == sp, bs)

                    reset_inbetween && Oscar.reset_global_serializer_state()
                    bs = load_glngraph_deform_bases(db, sp)
                    @test all(b -> b isa GlnGraphDeformBasis, bs)
                    @test length(bs) == 5
                    @test map(b -> b.degs, bs) == [0:0, 1:1, 2:2, 3:3, 4:4]
                    reset_inbetween || @test all(b -> b.sp == sp, bs)
                    reset_inbetween || @test collect(bs[3]) == collect(b_new)
                end

                @testset "loading PBW deformations; reset_inbetween=$(reset_inbetween)" begin
                    reset_inbetween && Oscar.reset_global_serializer_state()
                    ms = load_pbwdeformations(db, sp, 2:2)
                    @test ms isa Vector{<:DeformationMap}
                    reset_inbetween || @test all(m -> base_ring(m) == sp, ms)
                    ms_new = all_pbwdeformations(sp, GlnGraphDeformBasis(sp, 2:2))
                    @test length(ms) == length(ms_new)
                    reset_inbetween || @test ms == ms_new

                    reset_inbetween && Oscar.reset_global_serializer_state()
                    mss = load_pbwdeformations(db, sp, [0:3, 1:1])
                    @test all(ms -> ms isa Vector{<:DeformationMap}, mss)
                    reset_inbetween || @test all(ms -> all(m -> base_ring(m) == sp, ms), mss)
                    @test length(mss) == 2

                    reset_inbetween && Oscar.reset_global_serializer_state()
                    mss = load_pbwdeformations(db, sp)
                    @test all(ms -> ms isa Vector{<:DeformationMap}, mss)
                    reset_inbetween || @test all(ms -> all(m -> base_ring(m) == sp, ms), mss)
                    @test length(mss) == 5
                    reset_inbetween || @test mss[3] == ms_new
                end

            end
        end
    end
end
