"""
predictive_scores_tests.jl

測試預測分數評估功能
"""

using Test
using JSON2MAF

@testset "預測分數評估測試" begin

    # 創建預設配置
    config = FilterConfig()

    @testset "分數提取測試" begin
        # 測試完整的分數
        variant_full = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            0.85,  # primate_ai_3d
            0.80,  # primate_ai
            0.98,  # dann_score
            0.80,  # revel_score
            String[]
        )

        @test get_primate_ai_score(variant_full) == 0.85  # 優先使用 3D
        @test get_dann_score(variant_full) == 0.98
        @test get_revel_score(variant_full) == 0.80
        @test is_in_cosmic(variant_full) == false

        # 測試缺失分數
        variant_empty = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )

        @test get_primate_ai_score(variant_empty) === nothing
        @test get_dann_score(variant_empty) === nothing
        @test get_revel_score(variant_empty) === nothing
        @test is_in_cosmic(variant_empty) == false
    end

    @testset "PrimateAI 版本優先級" begin
        # 當兩個都存在時，優先使用 PrimateAI-3D
        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            0.85,  # primate_ai_3d
            0.70,  # primate_ai (舊版，較低)
            nothing, nothing,
            String[]
        )

        @test get_primate_ai_score(variant) == 0.85

        # 只有舊版時使用舊版
        variant_old = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing,  # primate_ai_3d 缺失
            0.70,     # primate_ai
            nothing, nothing,
            String[]
        )

        @test get_primate_ai_score(variant_old) == 0.70
    end

    @testset "COSMIC 檢測" begin
        cosmic_entry = CosmicEntry("COSM123", "TP53", "missense", 150)

        # 有 COSMIC 條目
        variant_with_cosmic = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], [cosmic_entry],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )

        @test is_in_cosmic(variant_with_cosmic) == true

        # 沒有 COSMIC 條目
        variant_no_cosmic = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )

        @test is_in_cosmic(variant_no_cosmic) == false
    end

    @testset "單一 PrimateAI-3D 支持" begin
        # PrimateAI-3D ≥ 0.8，單獨即可建議為 likely pathogenic
        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            0.85,  # primate_ai_3d ≥ 0.8
            nothing, nothing, nothing,
            String[]
        )

        assessment = assess_predictive_scores(variant, config)
        @test assessment.suggests_pathogenic == true
        @test haskey(assessment.contributing_scores, "PrimateAI-3D")
        @test assessment.contributing_scores["PrimateAI-3D"] == 0.85
        @test assessment.confidence >= 0.7  # PrimateAI-3D 基礎信心
    end

    @testset "2+ 分數支持" begin
        # REVEL + DANN 都達到閾值
        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing,  # 沒有 PrimateAI-3D
            nothing,
            0.98,     # DANN ≥ 0.96
            0.80,     # REVEL ≥ 0.75
            String[]
        )

        assessment = assess_predictive_scores(variant, config)
        @test assessment.suggests_pathogenic == true
        @test haskey(assessment.contributing_scores, "REVEL")
        @test haskey(assessment.contributing_scores, "DANN")
        @test length(assessment.contributing_scores) == 2
    end

    @testset "REVEL + COSMIC 支持" begin
        cosmic_entry = CosmicEntry("COSM123", "TP53", "missense", 150)

        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], [cosmic_entry],
            PopulationFrequency[],
            nothing,
            nothing,
            nothing,  # 沒有 DANN
            0.80,     # REVEL ≥ 0.75
            String[]
        )

        assessment = assess_predictive_scores(variant, config)
        @test assessment.suggests_pathogenic == true
        @test haskey(assessment.contributing_scores, "REVEL")
        @test haskey(assessment.contributing_scores, "COSMIC")
        @test assessment.contributing_scores["COSMIC"] == 1.0
    end

    @testset "單一分數不足" begin
        # 只有一個分數達標（非 PrimateAI-3D）
        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing,
            nothing,
            nothing,
            0.80,  # 只有 REVEL
            String[]
        )

        assessment = assess_predictive_scores(variant, config)
        @test assessment.suggests_pathogenic == false
        @test haskey(assessment.contributing_scores, "REVEL")
        @test length(assessment.contributing_scores) == 1
    end

    @testset "所有分數均低於閾值" begin
        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            0.70,  # PrimateAI-3D < 0.8
            nothing,
            0.90,  # DANN < 0.96
            0.60,  # REVEL < 0.75
            String[]
        )

        assessment = assess_predictive_scores(variant, config)
        @test assessment.suggests_pathogenic == false
        @test isempty(assessment.contributing_scores)
        @test assessment.confidence == 0.0
    end

    @testset "自訂配置閾值" begin
        # 使用較低的 REVEL 閾值 (0.5)
        config_low = create_filter_config(min_revel_score = 0.5)

        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing,
            0.98,  # DANN
            0.60,  # REVEL = 0.6 (≥ 0.5 但 < 0.75)
            String[]
        )

        # 使用預設配置 (0.75)，不支持
        assessment_default = assess_predictive_scores(variant, config)
        @test assessment_default.suggests_pathogenic == false

        # 使用低閾值配置 (0.5)，支持
        assessment_low = assess_predictive_scores(variant, config_low)
        @test assessment_low.suggests_pathogenic == true
        @test haskey(assessment_low.contributing_scores, "REVEL")
        @test haskey(assessment_low.contributing_scores, "DANN")
    end

    @testset "信心分數計算" begin
        # 沒有支持分數
        variant_none = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        assessment = assess_predictive_scores(variant_none, config)
        @test assessment.confidence == 0.0

        # 只有 PrimateAI-3D (1個)
        variant_primate = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            0.85, nothing, nothing, nothing,
            String[]
        )
        assessment = assess_predictive_scores(variant_primate, config)
        @test assessment.confidence == 0.7

        # PrimateAI-3D + REVEL (2個)
        variant_two = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            0.85, nothing, nothing, 0.80,
            String[]
        )
        assessment = assess_predictive_scores(variant_two, config)
        @test assessment.confidence ≈ 0.8

        # 所有分數都達標 (4個)
        cosmic_entry = CosmicEntry("COSM123", "TP53", "missense", 150)
        variant_all = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], [cosmic_entry],
            PopulationFrequency[],
            0.85, nothing, 0.98, 0.80,
            String[]
        )
        assessment = assess_predictive_scores(variant_all, config)
        @test assessment.confidence == 1.0  # 最高
    end

    @testset "複雜組合情境" begin
        cosmic_entry = CosmicEntry("COSM123", "TP53", "missense", 150)

        # PrimateAI-3D + REVEL + DANN + COSMIC
        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], [cosmic_entry],
            PopulationFrequency[],
            0.85,  # PrimateAI-3D
            nothing,
            0.98,  # DANN
            0.80,  # REVEL
            String[]
        )

        assessment = assess_predictive_scores(variant, config)
        @test assessment.suggests_pathogenic == true
        @test length(assessment.contributing_scores) == 4
        @test haskey(assessment.contributing_scores, "PrimateAI-3D")
        @test haskey(assessment.contributing_scores, "REVEL")
        @test haskey(assessment.contributing_scores, "DANN")
        @test haskey(assessment.contributing_scores, "COSMIC")
        @test assessment.confidence == 1.0
    end
end
