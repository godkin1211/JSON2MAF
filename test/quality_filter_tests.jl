"""
quality_filter_tests.jl

測試品質與族群頻率過濾功能
"""

using Test
using JSON2MAF

@testset "品質過濾器測試" begin
    # 創建預設配置
    config = FilterConfig()

    @testset "測序深度過濾" begin
        # 測試通過 - 深度剛好等於閾值
        variant_pass = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            30,  # totalDepth = 30 (剛好等於閾值)
            [0.5],  # VAF = 0.5
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_pass, config)
        @test result.passes_quality == true
        @test result.failure_reason === nothing

        # 測試失敗 - 深度低於閾值
        variant_fail = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            29,  # totalDepth = 29 (低於閾值)
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

        # 測試通過 - 深度高於閾值
        variant_high = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
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

    @testset "變異頻率(VAF)過濾" begin
        # 測試通過 - VAF 剛好等於閾值
        variant_pass = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100,
            [0.03],  # VAF = 0.03 (剛好等於閾值)
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_pass, config)
        @test result.passes_quality == true

        # 測試失敗 - VAF 低於閾值
        variant_fail = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100,
            [0.029],  # VAF = 0.029 (低於閾值)
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_fail, config)
        @test result.passes_quality == false
        @test occursin("Low variant frequency", result.failure_reason)
        @test occursin("0.029", result.failure_reason)

        # 測試通過 - VAF 高於閾值
        variant_high = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
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

    @testset "gnomAD-exome 東亞頻率過濾" begin
        # 測試通過 - easAf 低於閾值
        pop_freq_low = PopulationFrequency("gnomad-exome", 0.05, 0.005, nothing, nothing, nothing)
        variant_pass = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            [pop_freq_low],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_pass, config)
        @test result.passes_quality == true
        @test result.eas_allele_frequency == 0.005

        # 測試失敗 - easAf 高於閾值
        pop_freq_high = PopulationFrequency("gnomad-exome", 0.05, 0.02, nothing, nothing, nothing)
        variant_fail = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
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

        # 測試邊界值 - easAf 剛好等於閾值
        pop_freq_equal = PopulationFrequency("gnomad-exome", 0.05, 0.01, nothing, nothing, nothing)
        variant_equal = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            [pop_freq_equal],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_equal, config)
        @test result.passes_quality == true  # <= 閾值應該通過
    end

    @testset "1000 Genomes (oneKg) 東亞頻率過濾" begin
        # 測試通過 - easAf 低於閾值
        pop_freq_low = PopulationFrequency("oneKg", 0.05, 0.008, nothing, nothing, nothing)
        variant_pass = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            [pop_freq_low],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_pass, config)
        @test result.passes_quality == true
        @test result.eas_allele_frequency == 0.008

        # 測試失敗 - easAf 高於閾值
        pop_freq_high = PopulationFrequency("oneKg", 0.05, 0.015, nothing, nothing, nothing)
        variant_fail = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
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

    @testset "族群頻率優先順序" begin
        # gnomad-exome 優先於 oneKg
        pop_freq_gnomad = PopulationFrequency("gnomad-exome", 0.05, 0.005, nothing, nothing, nothing)
        pop_freq_onekg = PopulationFrequency("oneKg", 0.05, 0.02, nothing, nothing, nothing)

        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            [pop_freq_gnomad, pop_freq_onekg],  # gnomad-exome 在前
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant, config)
        # 應該使用 gnomad-exome 的值 (0.005)，因此通過
        @test result.passes_quality == true
        @test result.eas_allele_frequency == 0.005
    end

    @testset "缺失欄位處理" begin
        # 缺失 totalDepth
        variant_no_depth = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            nothing,  # 缺失 totalDepth
            [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_no_depth, config)
        @test result.passes_quality == false
        @test occursin("Missing totalDepth", result.failure_reason)

        # 缺失 variant_frequencies
        variant_no_vaf = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100,
            nothing,  # 缺失 VAF
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_no_vaf, config)
        @test result.passes_quality == false
        @test occursin("Missing variant frequency", result.failure_reason)

        # 缺失族群頻率 - 應該通過（保守策略）
        variant_no_pop = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],  # 沒有族群頻率資料
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_no_pop, config)
        @test result.passes_quality == true
        @test result.eas_allele_frequency === nothing
    end

    @testset "自訂配置參數" begin
        # 使用更嚴格的配置
        strict_config = create_filter_config(
            min_total_depth = 50,
            min_variant_frequency = 0.1,
            max_eas_af = 0.005
        )

        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            40,  # 不符合 strict_config 的深度要求
            [0.08],  # 不符合 strict_config 的 VAF 要求
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant, strict_config)
        @test result.passes_quality == false
        @test occursin("40", result.failure_reason)
    end

    @testset "多重條件組合" begin
        # 同時滿足所有條件
        pop_freq = PopulationFrequency("gnomad-exome", 0.05, 0.002, nothing, nothing, nothing)
        variant_all_pass = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            150,  # 高深度
            [0.45],  # 高 VAF
            TranscriptAnnotation[], ClinVarEntry[], CosmicEntry[],
            [pop_freq],  # 低族群頻率
            nothing, nothing, nothing, nothing,
            String[]
        )
        result = apply_quality_filters(variant_all_pass, config)
        @test result.passes_quality == true
        @test result.depth == 150
        @test result.variant_frequency == 0.45
        @test result.eas_allele_frequency == 0.002

        # 深度通過但 VAF 不通過
        variant_mixed = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            150,
            [0.01],  # VAF 過低
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
