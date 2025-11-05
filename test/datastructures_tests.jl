"""
datastructures_tests.jl

測試核心資料結構和配置模組
"""

using Test
using JSON2MAF

@testset "FilterConfig 測試" begin
    @testset "預設配置" begin
        config = FilterConfig()

        @test config.min_total_depth == 30
        @test config.min_variant_frequency == 0.03
        @test config.max_eas_af == 0.01
        @test config.min_revel_score == 0.75
        @test config.min_primate_ai_score == 0.8
        @test config.min_dann_score == 0.96
    end

    @testset "自訂配置" begin
        config = create_filter_config(
            min_total_depth = 50,
            max_eas_af = 0.005,
            min_revel_score = 0.5
        )

        @test config.min_total_depth == 50
        @test config.max_eas_af == 0.005
        @test config.min_revel_score == 0.5
        # 其他參數應該使用預設值
        @test config.min_variant_frequency == 0.03
        @test config.min_primate_ai_score == 0.8
    end

    @testset "配置驗證 - 合法值" begin
        # 這些應該不會拋出錯誤
        @test_nowarn create_filter_config(min_total_depth = 1)
        @test_nowarn create_filter_config(min_variant_frequency = 0.0)
        @test_nowarn create_filter_config(min_variant_frequency = 1.0)
        @test_nowarn create_filter_config(max_eas_af = 0.0)
        @test_nowarn create_filter_config(max_eas_af = 1.0)
        @test_nowarn create_filter_config(min_revel_score = 0.5)
        @test_nowarn create_filter_config(min_revel_score = 0.75)
    end

    @testset "配置驗證 - 非法值" begin
        # 深度不能小於 1
        @test_throws ErrorException create_filter_config(min_total_depth = 0)
        @test_throws ErrorException create_filter_config(min_total_depth = -1)

        # 頻率必須在 0-1 之間
        @test_throws ErrorException create_filter_config(min_variant_frequency = -0.1)
        @test_throws ErrorException create_filter_config(min_variant_frequency = 1.1)
        @test_throws ErrorException create_filter_config(max_eas_af = -0.1)
        @test_throws ErrorException create_filter_config(max_eas_af = 1.1)

        # 預測分數必須在 0-1 之間
        @test_throws ErrorException create_filter_config(min_revel_score = -0.1)
        @test_throws ErrorException create_filter_config(min_revel_score = 1.1)
        @test_throws ErrorException create_filter_config(min_primate_ai_score = 1.5)
        @test_throws ErrorException create_filter_config(min_dann_score = -0.5)
    end

    @testset "配置顯示" begin
        config = FilterConfig()
        # 測試 display_config 不會拋出錯誤
        io = IOBuffer()
        @test_nowarn display_config(config, io=io)

        # 檢查輸出包含關鍵資訊
        output = String(take!(io))
        @test occursin("JSON2MAF Filter Configuration", output)
        @test occursin("品質過濾參數", output)
        @test occursin("族群頻率過濾參數", output)
        @test occursin("預測分數閾值", output)
    end
end

@testset "資料結構初始化測試" begin
    @testset "ClinVarEntry" begin
        entry = ClinVarEntry(
            "RCV000123456",
            "12345",
            "Pathogenic",
            "criteria provided, multiple submitters, no conflicts",
            ["Hereditary cancer-predisposing syndrome"],
            "2021-10-12"
        )

        @test entry.id == "RCV000123456"
        @test entry.clinical_significance == "Pathogenic"
        @test length(entry.diseases) == 1
    end

    @testset "TranscriptAnnotation" begin
        transcript = TranscriptAnnotation(
            "NM_007294.3",
            "BRCA1",
            1100,
            ["missense_variant"],
            "V/M",
            1234,
            4321,
            988,
            "c.1234G>A",
            "p.Val988Met"
        )

        @test transcript.gene_symbol == "BRCA1"
        @test "missense_variant" in transcript.consequence
        @test transcript.hgvsc == "c.1234G>A"
    end

    @testset "VariantPosition 基本結構" begin
        variant = VariantPosition(
            "7",                          # chromosome
            140453136,                    # start
            140453136,                    # end_pos
            "G",                          # reference_allele
            "A",                          # alternate_allele
            "SNV",                        # variant_type
            100,                          # total_depth
            [0.45],                       # variant_frequencies
            TranscriptAnnotation[],       # transcripts
            ClinVarEntry[],              # clinvar
            CosmicEntry[],               # cosmic
            PopulationFrequency[],       # population_frequencies
            0.85,                         # primate_ai_3d
            0.82,                         # primate_ai
            0.97,                         # dann_score
            0.78,                         # revel_score
            String[]                      # dbsnp_ids
        )

        @test variant.chromosome == "7"
        @test variant.start == 140453136
        @test variant.variant_type == "SNV"
        @test variant.total_depth == 100
        @test variant.primate_ai_3d == 0.85
    end

    @testset "MAFRecord 基本結構" begin
        record = MAFRecord(
            hugo_symbol = "BRCA1",
            chromosome = "17",
            start_position = 41234567,
            end_position = 41234567,
            variant_classification = "Missense_Mutation",
            variant_type = "SNP",
            reference_allele = "G",
            tumor_seq_allele2 = "A"
        )

        @test record.hugo_symbol == "BRCA1"
        @test record.chromosome == "17"
        @test record.variant_classification == "Missense_Mutation"
    end

    @testset "FilterDecision" begin
        decision = FilterDecision(
            true,
            "Pathogenic",
            "ClinVar",
            "ClinVar Pathogenic with review status: criteria provided, multiple submitters"
        )

        @test decision.should_include == true
        @test decision.pathogenicity_class == "Pathogenic"
        @test decision.primary_evidence == "ClinVar"
    end
end

@testset "型別檢查" begin
    @testset "FilterConfig 型別" begin
        config = FilterConfig()
        @test config isa FilterConfig
        @test typeof(config.min_total_depth) == Int
        @test typeof(config.min_variant_frequency) == Float64
    end

    @testset "過濾結果型別" begin
        qf = QualityFilterResult(true, nothing, 50, 0.45, 0.005)
        @test qf isa QualityFilterResult
        @test qf.passes_quality == true

        ca = ClinVarAssessment(true, false, nothing, "high", "Pathogenic")
        @test ca isa ClinVarAssessment
        @test ca.is_pathogenic == true

        pa = PredictiveAssessment(true, Dict("REVEL" => 0.8), 0.8, 1, false)
        @test pa isa PredictiveAssessment
        @test pa.support_count == 1
    end
end
