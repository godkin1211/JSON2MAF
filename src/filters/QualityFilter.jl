"""
QualityFilter.jl

品質與族群頻率預過濾模組
在進行 ClinVar 評估前先過濾低品質與常見變異
"""

module QualityFilter

using ..DataStructures

export apply_quality_filters, check_sequencing_quality, check_population_frequency

"""
    apply_quality_filters(variant::VariantPosition, config::FilterConfig) -> QualityFilterResult

對單一變異位點應用品質過濾規則

# 過濾規則
1. 測序深度 (totalDepth) >= min_total_depth (預設 30)
2. 變異頻率 (variantFrequencies) >= min_variant_frequency (預設 0.03)
3. 東亞族群頻率 (easAf) <= max_eas_af (預設 0.01)

# 參數
- `variant`: 變異位點資料
- `config`: 過濾配置參數

# 返回
- `QualityFilterResult`: 包含是否通過及失敗原因

# 範例
```julia
result = apply_quality_filters(variant, config)
if result.passes_quality
    println("通過品質過濾")
else
    println("未通過: \$(result.failure_reason)")
end
```
"""
function apply_quality_filters(variant::VariantPosition, config::FilterConfig)::QualityFilterResult
    # 檢查 VCF filters 欄位：只接受 ["PASS"]
    if !(length(variant.filters) == 1 && variant.filters[1] == "PASS")
        filters_str = join(variant.filters, ", ")
        return QualityFilterResult(
            false,
            "Failed VCF filters: [$filters_str]",
            variant.total_depth,
            get_variant_frequency(variant),
            nothing
        )
    end

    # 檢查測序品質
    seq_pass, seq_reason = check_sequencing_quality(variant, config)
    if !seq_pass
        return QualityFilterResult(
            false,
            seq_reason,
            variant.total_depth,
            get_variant_frequency(variant),
            nothing
        )
    end

    # 檢查族群頻率
    pop_pass, pop_reason, eas_af = check_population_frequency(variant, config)
    if !pop_pass
        return QualityFilterResult(
            false,
            pop_reason,
            variant.total_depth,
            get_variant_frequency(variant),
            eas_af
        )
    end

    # 全部通過
    return QualityFilterResult(
        true,
        nothing,
        variant.total_depth,
        get_variant_frequency(variant),
        eas_af
    )
end

"""
    check_sequencing_quality(variant::VariantPosition, config::FilterConfig) -> Tuple{Bool, String}

檢查測序品質是否符合標準

# 檢查項目
- totalDepth >= min_total_depth
- variantFrequencies >= min_variant_frequency

# 返回
- (通過與否, 失敗原因)
"""
function check_sequencing_quality(variant::VariantPosition, config::FilterConfig)::Tuple{Bool, String}
    # 檢查測序深度
    if variant.total_depth === nothing
        return (false, "Missing totalDepth")
    end

    if variant.total_depth < config.min_total_depth
        return (false, "Low sequencing depth ($(variant.total_depth) < $(config.min_total_depth))")
    end

    # 檢查變異頻率
    vaf = get_variant_frequency(variant)
    if vaf === nothing
        return (false, "Missing variant frequency")
    end

    if vaf < config.min_variant_frequency
        return (false, "Low variant frequency ($(round(vaf, digits=4)) < $(config.min_variant_frequency))")
    end

    return (true, "")
end

"""
    check_population_frequency(variant::VariantPosition, config::FilterConfig) -> Tuple{Bool, String, Union{Float64, Nothing}}

檢查東亞族群頻率是否低於閾值

優先順序:
1. 檢查 gnomad-exome 的 easAf
2. 檢查 oneKg 的 easAf

# 返回
- (通過與否, 失敗原因, 東亞AF值)
"""
function check_population_frequency(variant::VariantPosition, config::FilterConfig)::Tuple{Bool, String, Union{Float64, Nothing}}
    # 嘗試從 gnomad-exome 提取 easAf
    gnomad_exome_af = extract_gnomad_exome_eas_af(variant)
    if gnomad_exome_af !== nothing
        if gnomad_exome_af > config.max_eas_af
            return (false, "High East Asian AF in gnomAD-exome ($(round(gnomad_exome_af, digits=4)) > $(config.max_eas_af))", gnomad_exome_af)
        end
        return (true, "", gnomad_exome_af)
    end

    # 嘗試從 oneKg 提取 easAf
    onekg_af = extract_onekg_eas_af(variant)
    if onekg_af !== nothing
        if onekg_af > config.max_eas_af
            return (false, "High East Asian AF in 1000G ($(round(onekg_af, digits=4)) > $(config.max_eas_af))", onekg_af)
        end
        return (true, "", onekg_af)
    end

    # 沒有族群頻率資料，視為通過（保守策略）
    return (true, "", nothing)
end

"""
    extract_gnomad_exome_eas_af(variant::VariantPosition) -> Union{Float64, Nothing}

從 gnomad-exome 族群頻率中提取東亞族群 AF
"""
function extract_gnomad_exome_eas_af(variant::VariantPosition)::Union{Float64, Nothing}
    for pf in variant.population_frequencies
        if pf.source == "gnomad-exome" && pf.eas_af !== nothing
            return pf.eas_af
        end
    end
    return nothing
end

"""
    extract_onekg_eas_af(variant::VariantPosition) -> Union{Float64, Nothing}

從 1000 Genomes 族群頻率中提取東亞族群 AF
"""
function extract_onekg_eas_af(variant::VariantPosition)::Union{Float64, Nothing}
    for pf in variant.population_frequencies
        if pf.source == "oneKg" && pf.eas_af !== nothing
            return pf.eas_af
        end
    end
    return nothing
end

"""
    get_variant_frequency(variant::VariantPosition) -> Union{Float64, Nothing}

提取變異頻率（取第一個值）
"""
function get_variant_frequency(variant::VariantPosition)::Union{Float64, Nothing}
    if variant.variant_frequencies !== nothing && length(variant.variant_frequencies) > 0
        return variant.variant_frequencies[1]
    end
    return nothing
end

end # module QualityFilter
