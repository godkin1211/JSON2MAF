"""
quality_filter_tests.jl

Test quality and population frequency filtering functionality
"""

using Test
using JSON2MAF

@testset "Quality Filter Tests" begin
    # Create default configuration
    config = FilterConfig()

    @testset "Sequencing Depth Filtering" begin
        # Test pass - depth exactly equals threshold
        variant_pass = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            30,  # totalDepth = 30 (exactly equals threshold)
            [0.5],  # VAF = 0.5
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_pass, config)
        @test result.passes_quality == true
        @test result.failure_reason === nothing

        # Test fail - depth below threshold
        variant_fail = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            29,  # totalDepth = 29 (below threshold)
            [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_fail, config)
        @test result.passes_quality == false
        @test occursin("Low sequencing depth", result.failure_reason)
        @test occursin("29", result.failure_reason)

        # Test pass - depth above threshold
        variant_high = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100,
            [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_high, config)
        @test result.passes_quality == true
    end

    @testset "Variant Frequency (VAF) Filtering" begin
        # Test pass - VAF exactly equals threshold
        variant_pass = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100,
            [0.03],  # VAF = 0.03 (exactly equals threshold)
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_pass, config)
        @test result.passes_quality == true

        # Test fail - VAF below threshold
        variant_fail = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100,
            [0.029],  # VAF = 0.029 (below threshold)
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_fail, config)
        @test result.passes_quality == false
        @test occursin("Low variant frequency", result.failure_reason)
        @test occursin("0.029", result.failure_reason)

        # Test pass - VAF above threshold
        variant_high = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100,
            [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_high, config)
        @test result.passes_quality == true
    end

    @testset "gnomAD-exome East Asian Frequency Filtering" begin
        # Test pass - easAf below threshold
        pop_freq_low = PopulationFrequency("gnomad-exome", 0.05, 0.005, nothing, nothing, nothing)
        variant_pass = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            [pop_freq_low],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_pass, config)
        @test result.passes_quality == true
        @test result.eas_allele_frequency == 0.005

        # Test fail - easAf above threshold
        pop_freq_high = PopulationFrequency("gnomad-exome", 0.05, 0.02, nothing, nothing, nothing)
        variant_fail = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            [pop_freq_high],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_fail, config)
        @test result.passes_quality == false
        @test occursin("High East Asian AF", result.failure_reason)
        @test occursin("gnomAD-exome", result.failure_reason)
        @test result.eas_allele_frequency == 0.02

        # Test boundary value - easAf exactly equals threshold
        pop_freq_equal = PopulationFrequency("gnomad-exome", 0.05, 0.01, nothing, nothing, nothing)
        variant_equal = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            [pop_freq_equal],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_equal, config)
        @test result.passes_quality == true  # <= threshold should pass
    end

    @testset "1000 Genomes (oneKg) East Asian Frequency Filtering" begin
        # Test pass - easAf below threshold
        pop_freq_low = PopulationFrequency("oneKg", 0.05, 0.008, nothing, nothing, nothing)
        variant_pass = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            [pop_freq_low],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_pass, config)
        @test result.passes_quality == true
        @test result.eas_allele_frequency == 0.008

        # Test fail - easAf above threshold
        pop_freq_high = PopulationFrequency("oneKg", 0.05, 0.015, nothing, nothing, nothing)
        variant_fail = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            [pop_freq_high],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_fail, config)
        @test result.passes_quality == false
        @test occursin("High East Asian AF", result.failure_reason)
        @test occursin("1000G", result.failure_reason)
    end

    @testset "Population Frequency Priority Order" begin
        # gnomad-exome takes priority over oneKg
        pop_freq_gnomad = PopulationFrequency("gnomad-exome", 0.05, 0.005, nothing, nothing, nothing)
        pop_freq_onekg = PopulationFrequency("oneKg", 0.05, 0.02, nothing, nothing, nothing)

        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            [pop_freq_gnomad, pop_freq_onekg],  # gnomad-exome first
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant, config)
        # Should use gnomad-exome value (0.005), therefore passes
        @test result.passes_quality == true
        @test result.eas_allele_frequency == 0.005
    end

    @testset "Missing Field Handling" begin
        # Missing totalDepth
        variant_no_depth = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            nothing,  # Missing totalDepth
            [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_no_depth, config)
        @test result.passes_quality == false
        @test occursin("Missing totalDepth", result.failure_reason)

        # Missing variant_frequencies
        variant_no_vaf = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100,
            nothing,  # Missing VAF
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_no_vaf, config)
        @test result.passes_quality == false
        @test occursin("Missing variant frequency", result.failure_reason)

        # Missing population frequency - should pass (conservative strategy)
        variant_no_pop = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],  # No population frequency data
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_no_pop, config)
        @test result.passes_quality == true
        @test result.eas_allele_frequency === nothing
    end

    @testset "Custom Configuration Parameters" begin
        # Use stricter configuration
        strict_config = create_filter_config(
            min_total_depth = 50,
            min_variant_frequency = 0.1,
            max_eas_af = 0.005
        )

        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            40,  # Does not meet strict_config depth requirement
            [0.08],  # Does not meet strict_config VAF requirement
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant, strict_config)
        @test result.passes_quality == false
        @test occursin("40", result.failure_reason)
    end

    @testset "Multiple Condition Combinations" begin
        # Satisfies all conditions
        pop_freq = PopulationFrequency("gnomad-exome", 0.05, 0.002, nothing, nothing, nothing)
        variant_all_pass = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            150,  # High depth
            [0.45],  # High VAF
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            [pop_freq],  # Low population frequency
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_all_pass, config)
        @test result.passes_quality == true
        @test result.depth == 150
        @test result.variant_frequency == 0.45
        @test result.eas_allele_frequency == 0.002

        # Depth passes but VAF does not pass
        variant_mixed = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV", ["PASS"],
            150,
            [0.01],  # VAF too low
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            [pop_freq],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_mixed, config)
        @test result.passes_quality == false
        @test occursin("variant frequency", result.failure_reason)
    end
end
