"""
QualityFilter.jl

Quality and population frequency prefiltering module
Filters low-quality and common variants before ClinVar assessment
"""

module QualityFilter

using ..DataStructures

export apply_quality_filters, check_sequencing_quality, check_population_frequency

"""
    apply_quality_filters(variant::VariantPosition, config::FilterConfig) -> QualityFilterResult

Apply quality filtering rules to a single variant position

# Filtering rules
1. Sequencing depth (totalDepth) >= min_total_depth (default 30)
2. Variant frequency (variantFrequencies) >= min_variant_frequency (default 0.03)
3. East Asian population frequency (easAf) <= max_eas_af (default 0.01)

# Parameters
- `variant`: Variant position data
- `config`: Filter configuration parameters

# Returns
- `QualityFilterResult`: Contains whether passed and failure reason

# Example
```julia
result = apply_quality_filters(variant, config)
if result.passes_quality
    println("Passed quality filtering")
else
    println("Failed: \$(result.failure_reason)")
end
```
"""
function apply_quality_filters(variant::VariantPosition, config::FilterConfig)::QualityFilterResult
    # Check VCF filters field: only accept ["PASS"]
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

    # Check sequencing quality
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

    # Check population frequency
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

    # All passed
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

Check if sequencing quality meets standards

# Check items
- totalDepth >= min_total_depth
- variantFrequencies >= min_variant_frequency

# Returns
- (Pass or fail, Failure reason)
"""
function check_sequencing_quality(variant::VariantPosition, config::FilterConfig)::Tuple{Bool, String}
    # Check sequencing depth
    if variant.total_depth === nothing
        return (false, "Missing totalDepth")
    end

    if variant.total_depth < config.min_total_depth
        return (false, "Low sequencing depth ($(variant.total_depth) < $(config.min_total_depth))")
    end

    # Check variant frequency
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

Check if East Asian population frequency is below threshold

Priority:
1. Check gnomad-exome easAf
2. Check oneKg easAf

# Returns
- (Pass or fail, Failure reason, East Asian AF value)
"""
function check_population_frequency(variant::VariantPosition, config::FilterConfig)::Tuple{Bool, String, Union{Float64, Nothing}}
    # Try to extract easAf from gnomad-exome
    gnomad_exome_af = extract_gnomad_exome_eas_af(variant)
    if gnomad_exome_af !== nothing
        if gnomad_exome_af > config.max_eas_af
            return (false, "High East Asian AF in gnomAD-exome ($(round(gnomad_exome_af, digits=4)) > $(config.max_eas_af))", gnomad_exome_af)
        end
        return (true, "", gnomad_exome_af)
    end

    # Try to extract easAf from oneKg
    onekg_af = extract_onekg_eas_af(variant)
    if onekg_af !== nothing
        if onekg_af > config.max_eas_af
            return (false, "High East Asian AF in 1000G ($(round(onekg_af, digits=4)) > $(config.max_eas_af))", onekg_af)
        end
        return (true, "", onekg_af)
    end

    # No population frequency data, consider as pass (conservative strategy)
    return (true, "", nothing)
end

"""
    extract_gnomad_exome_eas_af(variant::VariantPosition) -> Union{Float64, Nothing}

Extract East Asian population AF from gnomad-exome population frequency
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

Extract East Asian population AF from 1000 Genomes population frequency
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

Extract variant frequency (take first value)
"""
function get_variant_frequency(variant::VariantPosition)::Union{Float64, Nothing}
    if variant.variant_frequencies !== nothing && length(variant.variant_frequencies) > 0
        return variant.variant_frequencies[1]
    end
    return nothing
end

end # module QualityFilter
