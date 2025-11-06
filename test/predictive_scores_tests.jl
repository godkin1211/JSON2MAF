"""
predictive_scores_tests.jl

Test predictive score assessment functionality
"""

using Test
using JSON2MAF

@testset "Predictive Score Assessment Tests" begin

    # Create default configuration
    config = FilterConfig()

    @testset "Score Extraction Tests" begin
        # Test complete scores
        variant_full = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            0.85,  # primate_ai_3d
            0.80,  # primate_ai
            0.98,  # dann_score
            0.80,  # revel_score
            String[]
        )

        @test get_primate_ai_score(variant_full) == 0.85  # Prioritize 3D
        @test get_dann_score(variant_full) == 0.98
        @test get_revel_score(variant_full) == 0.80
        @test is_in_cosmic(variant_full) == false

        # Test missing scores
        variant_empty = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
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

    @testset "PrimateAI Version Priority" begin
        # When both exist, prioritize PrimateAI-3D
        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            0.85,  # primate_ai_3d
            0.70,  # primate_ai (old version, lower)
            nothing, nothing,
            String[]
        )

        @test get_primate_ai_score(variant) == 0.85

        # Use old version when only old version exists
        variant_old = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing,  # primate_ai_3d missing
            0.70,     # primate_ai
            nothing, nothing,
            String[]
        )

        @test get_primate_ai_score(variant_old) == 0.70
    end

    @testset "COSMIC Detection" begin
        cosmic_entry = CosmicEntry("COSM123", "TP53", "missense", 150)

        # Has COSMIC entry
        variant_with_cosmic = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], [cosmic_entry],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )

        @test is_in_cosmic(variant_with_cosmic) == true

        # No COSMIC entry
        variant_no_cosmic = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )

        @test is_in_cosmic(variant_no_cosmic) == false
    end

    @testset "Single PrimateAI-3D Support" begin
        # PrimateAI-3D ≥ 0.8, alone can suggest likely pathogenic
        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            0.85,  # primate_ai_3d >= 0.8
            nothing, nothing, nothing,
            String[]
        )

        assessment = assess_predictive_scores(variant, config)
        @test assessment.suggests_pathogenic == true
        @test haskey(assessment.contributing_scores, "PrimateAI-3D")
        @test assessment.contributing_scores["PrimateAI-3D"] == 0.85
        @test assessment.confidence >= 0.7  # PrimateAI-3D base confidence
    end

    @testset "2+ Score Support" begin
        # REVEL + DANN both meet threshold
        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing,  # No PrimateAI-3D
            nothing,
            0.98,     # DANN >= 0.96
            0.80,     # REVEL >= 0.75
            String[]
        )

        assessment = assess_predictive_scores(variant, config)
        @test assessment.suggests_pathogenic == true
        @test haskey(assessment.contributing_scores, "REVEL")
        @test haskey(assessment.contributing_scores, "DANN")
        @test length(assessment.contributing_scores) == 2
    end

    @testset "REVEL + COSMIC Support" begin
        cosmic_entry = CosmicEntry("COSM123", "TP53", "missense", 150)

        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], [cosmic_entry],
            PopulationFrequency[],
            nothing,
            nothing,
            nothing,  # No DANN
            0.80,     # REVEL >= 0.75
            String[]
        )

        assessment = assess_predictive_scores(variant, config)
        @test assessment.suggests_pathogenic == true
        @test haskey(assessment.contributing_scores, "REVEL")
        @test haskey(assessment.contributing_scores, "COSMIC")
        @test assessment.contributing_scores["COSMIC"] == 1.0
    end

    @testset "Single Score Insufficient" begin
        # Only one score meets threshold (not PrimateAI-3D)
        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing,
            nothing,
            nothing,
            0.80,  # Only REVEL
            String[]
        )

        assessment = assess_predictive_scores(variant, config)
        @test assessment.suggests_pathogenic == false
        @test haskey(assessment.contributing_scores, "REVEL")
        @test length(assessment.contributing_scores) == 1
    end

    @testset "All Scores Below Threshold" begin
        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
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

    @testset "Custom Configuration Thresholds" begin
        # Use lower REVEL threshold (0.5)
        config_low = create_filter_config(min_revel_score = 0.5)

        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing,
            0.98,  # DANN
            0.60,  # REVEL = 0.6 (>= 0.5 but < 0.75)
            String[]
        )

        # Use default config (0.75), not supported
        assessment_default = assess_predictive_scores(variant, config)
        @test assessment_default.suggests_pathogenic == false

        # Use low threshold config (0.5), supported
        assessment_low = assess_predictive_scores(variant, config_low)
        @test assessment_low.suggests_pathogenic == true
        @test haskey(assessment_low.contributing_scores, "REVEL")
        @test haskey(assessment_low.contributing_scores, "DANN")
    end

    @testset "Confidence Score Calculation" begin
        # No supporting scores
        variant_none = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        assessment = assess_predictive_scores(variant_none, config)
        @test assessment.confidence == 0.0

        # Only PrimateAI-3D (1 score)
        variant_primate = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            0.85, nothing, nothing, nothing,
            String[]
        )
        assessment = assess_predictive_scores(variant_primate, config)
        @test assessment.confidence == 0.7

        # PrimateAI-3D + REVEL (2 scores)
        variant_two = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            0.85, nothing, nothing, 0.80,
            String[]
        )
        assessment = assess_predictive_scores(variant_two, config)
        @test assessment.confidence ≈ 0.8

        # All scores meet threshold (4 scores)
        cosmic_entry = CosmicEntry("COSM123", "TP53", "missense", 150)
        variant_all = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], [cosmic_entry],
            PopulationFrequency[],
            0.85, nothing, 0.98, 0.80,
            String[]
        )
        assessment = assess_predictive_scores(variant_all, config)
        @test assessment.confidence == 1.0  # Highest
    end

    @testset "Complex Combination Scenarios" begin
        cosmic_entry = CosmicEntry("COSM123", "TP53", "missense", 150)

        # PrimateAI-3D + REVEL + DANN + COSMIC
        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
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
