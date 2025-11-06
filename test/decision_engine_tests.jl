"""
decision_engine_tests.jl

Test integrated filtering decision engine
"""

using Test
using JSON2MAF

@testset "Integrated Filter Decision Engine Tests" begin

    # Create default configuration
    config = FilterConfig()

    # Helper function: Create test variant
    function create_test_variant(;primate_ai_3d=nothing, dann=nothing, revel=nothing, cosmic=false)
        cosmic_entries = cosmic ? [CosmicEntry("COSM123", "TP53", "missense", 150)] : CosmicEntry[]
        return VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], cosmic_entries,
            PopulationFrequency[],
            primate_ai_3d, nothing, dann, revel,
            String[]
        )
    end

    @testset "Rule 1: ClinVar Pathogenic → Direct Inclusion" begin
        variant = create_test_variant()

        clinvar_entry = ClinVarEntry(
            "RCV12345", "12345",
            "Pathogenic",
            "reviewed by expert panel",
            ["Breast cancer"],
            "2024-01-01"
        )
        clinvar_assessment = assess_clinvar_pathogenicity([clinvar_entry])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == true
        @test decision.pathogenicity_class == "Pathogenic"
        @test decision.primary_evidence == "ClinVar"
        @test occursin("ClinVar: Pathogenic", decision.justification)
        @test occursin("expert panel", decision.justification)
    end

    @testset "Rule 2: ClinVar Likely Pathogenic → Direct Inclusion" begin
        variant = create_test_variant()

        clinvar_entry = ClinVarEntry(
            "RCV67890", "67890",
            "Likely pathogenic",
            "criteria provided, multiple submitters, no conflicts",
            ["Cancer syndrome"],
            "2024-01-01"
        )
        clinvar_assessment = assess_clinvar_pathogenicity([clinvar_entry])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == true
        @test decision.pathogenicity_class == "Likely pathogenic"
        @test decision.primary_evidence == "ClinVar"
        @test occursin("Likely pathogenic", decision.justification)
    end

    @testset "Rule 3: ClinVar Inconclusive + 2+ Predictive Scores → Include" begin
        # REVEL + DANN
        variant = create_test_variant(revel=0.80, dann=0.98)

        clinvar_assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == true
        @test decision.pathogenicity_class == "Likely pathogenic"
        @test decision.primary_evidence == "Predictive"
        @test occursin("2 predictive scores", decision.justification)
        @test occursin("DANN", decision.justification)
        @test occursin("REVEL", decision.justification)
    end

    @testset "Rule 3: PrimateAI-3D + REVEL → Include" begin
        variant = create_test_variant(primate_ai_3d=0.85, revel=0.80)

        clinvar_assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == true
        @test decision.pathogenicity_class == "Likely pathogenic"
        @test decision.primary_evidence == "Predictive"
        @test occursin("2 predictive scores", decision.justification)
    end

    @testset "Rule 3: REVEL + DANN + COSMIC → Include" begin
        variant = create_test_variant(revel=0.80, dann=0.98, cosmic=true)

        clinvar_assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == true
        @test decision.pathogenicity_class == "Likely pathogenic"
        @test decision.primary_evidence == "Predictive"
        @test occursin("3 predictive scores", decision.justification)
    end

    @testset "Rule 4: Only PrimateAI-3D Support → Include (Highest Priority)" begin
        variant = create_test_variant(primate_ai_3d=0.85)

        clinvar_assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == true
        @test decision.pathogenicity_class == "Likely pathogenic"
        @test decision.primary_evidence == "Predictive"
        @test occursin("PrimateAI-3D", decision.justification)
        @test occursin("highest priority", decision.justification)
    end

    @testset "Rule 5: Only REVEL Support → Exclude (Insufficient Evidence)" begin
        variant = create_test_variant(revel=0.80)

        clinvar_assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == false
        @test decision.pathogenicity_class == "Excluded"
        @test decision.primary_evidence == "Insufficient evidence"
        @test occursin("REVEL", decision.justification)
        @test occursin("insufficient", decision.justification)
    end

    @testset "Rule 5: Only DANN Support → Exclude" begin
        variant = create_test_variant(dann=0.98)

        clinvar_assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == false
        @test decision.pathogenicity_class == "Excluded"
        @test occursin("DANN", decision.justification)
    end

    @testset "Rule 5: Only COSMIC Support → Exclude" begin
        variant = create_test_variant(cosmic=true)

        clinvar_assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == false
        @test decision.pathogenicity_class == "Excluded"
        @test occursin("COSMIC", decision.justification)
    end

    @testset "Rule 5: No Support → Exclude" begin
        variant = create_test_variant()

        clinvar_assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == false
        @test decision.pathogenicity_class == "Excluded"
        @test occursin("No pathogenic ClinVar", decision.justification)
        @test occursin("No supporting predictive", decision.justification)
    end

    @testset "Helper Function Tests" begin
        # count_supporting_predictive_scores
        variant = create_test_variant(revel=0.80, dann=0.98)
        predictive_assessment = assess_predictive_scores(variant, config)

        count = count_supporting_predictive_scores(predictive_assessment)
        @test count == 2

        # has_primate_ai_support
        variant_with_primate = create_test_variant(primate_ai_3d=0.85)
        assessment_with = assess_predictive_scores(variant_with_primate, config)
        @test has_primate_ai_support(assessment_with) == true

        variant_without_primate = create_test_variant(revel=0.80)
        assessment_without = assess_predictive_scores(variant_without_primate, config)
        @test has_primate_ai_support(assessment_without) == false
    end

    @testset "ClinVar Takes Priority Over Predictive Scores" begin
        # Even with multiple predictive scores, ClinVar Pathogenic should take priority
        variant = create_test_variant(primate_ai_3d=0.85, revel=0.80, dann=0.98)

        clinvar_entry = ClinVarEntry(
            "RCV12345", "12345",
            "Pathogenic",
            "reviewed by expert panel",
            ["Cancer"],
            "2024-01-01"
        )
        clinvar_assessment = assess_clinvar_pathogenicity([clinvar_entry])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == true
        @test decision.pathogenicity_class == "Pathogenic"  # Not Likely pathogenic
        @test decision.primary_evidence == "ClinVar"
    end

    @testset "Benign ClinVar + Strong Predictive Scores → Exclude (ClinVar Priority)" begin
        # When ClinVar is Benign, exclude even if predictive scores support
        variant = create_test_variant(primate_ai_3d=0.85, revel=0.80)

        benign_entry = ClinVarEntry(
            "RCV99999", "99999",
            "Benign",
            "reviewed by expert panel",
            ["Disease"],
            "2024-01-01"
        )
        clinvar_assessment = assess_clinvar_pathogenicity([benign_entry])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        # ClinVar Benign will be filtered by assess_clinvar_pathogenicity
        # Therefore clinvar_assessment will not have pathogenic evidence
        # Should be included based on predictive scores
        @test decision.should_include == true
        @test decision.primary_evidence == "Predictive"
    end

    @testset "Complex Scenario: All Assessment Combinations" begin
        # All 4 scores met + ClinVar Pathogenic
        variant = create_test_variant(primate_ai_3d=0.85, revel=0.80, dann=0.98, cosmic=true)

        clinvar_entry = ClinVarEntry(
            "RCV12345", "12345",
            "Pathogenic",
            "practice guideline",
            ["Hereditary cancer"],
            "2024-01-01"
        )
        clinvar_assessment = assess_clinvar_pathogenicity([clinvar_entry])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == true
        @test decision.pathogenicity_class == "Pathogenic"
        @test decision.primary_evidence == "ClinVar"
        @test occursin("high", decision.justification)  # practice guideline = high confidence
    end
end
