"""
decision_engine_tests.jl

測試整合過濾決策引擎
"""

using Test
using JSON2MAF

@testset "整合過濾決策引擎測試" begin

    # 創建預設配置
    config = FilterConfig()

    # 輔助函數：創建測試用變異
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

    @testset "規則 1: ClinVar Pathogenic → 直接納入" begin
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

    @testset "規則 2: ClinVar Likely pathogenic → 直接納入" begin
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

    @testset "規則 3: ClinVar 無結論 + 2+ 預測分數 → 納入" begin
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

    @testset "規則 3: PrimateAI-3D + REVEL → 納入" begin
        variant = create_test_variant(primate_ai_3d=0.85, revel=0.80)

        clinvar_assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == true
        @test decision.pathogenicity_class == "Likely pathogenic"
        @test decision.primary_evidence == "Predictive"
        @test occursin("2 predictive scores", decision.justification)
    end

    @testset "規則 3: REVEL + DANN + COSMIC → 納入" begin
        variant = create_test_variant(revel=0.80, dann=0.98, cosmic=true)

        clinvar_assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == true
        @test decision.pathogenicity_class == "Likely pathogenic"
        @test decision.primary_evidence == "Predictive"
        @test occursin("3 predictive scores", decision.justification)
    end

    @testset "規則 4: 僅 PrimateAI-3D 支持 → 納入（最高優先級）" begin
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

    @testset "規則 5: 僅 REVEL 支持 → 排除（證據不足）" begin
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

    @testset "規則 5: 僅 DANN 支持 → 排除" begin
        variant = create_test_variant(dann=0.98)

        clinvar_assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == false
        @test decision.pathogenicity_class == "Excluded"
        @test occursin("DANN", decision.justification)
    end

    @testset "規則 5: 僅 COSMIC 支持 → 排除" begin
        variant = create_test_variant(cosmic=true)

        clinvar_assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == false
        @test decision.pathogenicity_class == "Excluded"
        @test occursin("COSMIC", decision.justification)
    end

    @testset "規則 5: 無任何支持 → 排除" begin
        variant = create_test_variant()

        clinvar_assessment = assess_clinvar_pathogenicity(ClinVarEntry[])
        predictive_assessment = assess_predictive_scores(variant, config)

        decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

        @test decision.should_include == false
        @test decision.pathogenicity_class == "Excluded"
        @test occursin("No pathogenic ClinVar", decision.justification)
        @test occursin("No supporting predictive", decision.justification)
    end

    @testset "輔助函數測試" begin
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

    @testset "ClinVar 優先於預測分數" begin
        # 即使有多個預測分數，ClinVar Pathogenic 仍應優先
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
        @test decision.pathogenicity_class == "Pathogenic"  # 不是 Likely pathogenic
        @test decision.primary_evidence == "ClinVar"
    end

    @testset "Benign ClinVar + 強預測分數 → 排除（ClinVar 優先）" begin
        # ClinVar 為 Benign 時，即使預測分數支持，仍排除
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

        # ClinVar Benign 會被 assess_clinvar_pathogenicity 過濾掉
        # 因此 clinvar_assessment 不會有 pathogenic 證據
        # 應該根據預測分數納入
        @test decision.should_include == true
        @test decision.primary_evidence == "Predictive"
    end

    @testset "複雜情境：所有評估組合" begin
        # 4個分數都達標 + ClinVar Pathogenic
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
