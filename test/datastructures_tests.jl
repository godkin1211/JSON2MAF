"""
datastructures_tests.jl

Test core data structures and configuration modules
"""

using Test
using JSON2MAF

@testset "FilterConfig Tests" begin
    @testset "Default Configuration" begin
        config = FilterConfig()

        @test config.min_total_depth == 30
        @test config.min_variant_frequency == 0.03
        @test config.max_eas_af == 0.01
        @test config.min_revel_score == 0.75
        @test config.min_primate_ai_score == 0.8
        @test config.min_dann_score == 0.96
    end

    @testset "Custom Configuration" begin
        config = create_filter_config(
            min_total_depth = 50,
            max_eas_af = 0.005,
            min_revel_score = 0.5
        )

        @test config.min_total_depth == 50
        @test config.max_eas_af == 0.005
        @test config.min_revel_score == 0.5
        # Other parameters should use default values
        @test config.min_variant_frequency == 0.03
        @test config.min_primate_ai_score == 0.8
    end

    @testset "Configuration Validation - Valid Values" begin
        # These should not throw errors
        @test_nowarn create_filter_config(min_total_depth = 1)
        @test_nowarn create_filter_config(min_variant_frequency = 0.0)
        @test_nowarn create_filter_config(min_variant_frequency = 1.0)
        @test_nowarn create_filter_config(max_eas_af = 0.0)
        @test_nowarn create_filter_config(max_eas_af = 1.0)
        @test_nowarn create_filter_config(min_revel_score = 0.5)
        @test_nowarn create_filter_config(min_revel_score = 0.75)
    end

    @testset "Configuration Validation - Invalid Values" begin
        # Depth cannot be less than 1
        @test_throws ErrorException create_filter_config(min_total_depth = 0)
        @test_throws ErrorException create_filter_config(min_total_depth = -1)

        # Frequency must be between 0-1
        @test_throws ErrorException create_filter_config(min_variant_frequency = -0.1)
        @test_throws ErrorException create_filter_config(min_variant_frequency = 1.1)
        @test_throws ErrorException create_filter_config(max_eas_af = -0.1)
        @test_throws ErrorException create_filter_config(max_eas_af = 1.1)

        # Predictive scores must be between 0-1
        @test_throws ErrorException create_filter_config(min_revel_score = -0.1)
        @test_throws ErrorException create_filter_config(min_revel_score = 1.1)
        @test_throws ErrorException create_filter_config(min_primate_ai_score = 1.5)
        @test_throws ErrorException create_filter_config(min_dann_score = -0.5)
    end

    @testset "Configuration Display" begin
        config = FilterConfig()
        # Test that display_config does not throw errors
        io = IOBuffer()
        @test_nowarn display_config(config, io=io)

        # Check that output contains key information
        output = String(take!(io))
        @test occursin("JSON2MAF Filter Configuration", output)
        @test occursin("Quality filtering parameters", output)
        @test occursin("Population frequency filtering parameters", output)
        @test occursin("Predictive score thresholds", output)
    end
end

@testset "Data Structure Initialization Tests" begin
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
            "p.Val988Met",
            true  # is_mane_select
        )

        @test transcript.gene_symbol == "BRCA1"
        @test "missense_variant" in transcript.consequence
        @test transcript.hgvsc == "c.1234G>A"
    end

    @testset "VariantPosition Basic Structure" begin
        variant = VariantPosition(
            "7",                          # chromosome
            140453136,                    # start
            140453136,                    # end_pos
            "G",                          # reference_allele
            "A",                          # alternate_allele
            "SNV",                        # variant_type
            ["PASS"],                     # filters
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

    @testset "MAFRecord Basic Structure" begin
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

@testset "Type Checking" begin
    @testset "FilterConfig Type" begin
        config = FilterConfig()
        @test config isa FilterConfig
        @test typeof(config.min_total_depth) == Int
        @test typeof(config.min_variant_frequency) == Float64
    end

    @testset "Filter Result Types" begin
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
