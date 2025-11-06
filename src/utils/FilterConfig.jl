"""
FilterConfig.jl

Filter configuration management module
Provides configuration creation, validation, and display functions
"""

module FilterConfigModule

using ..DataStructures

export create_filter_config, validate_config, display_config

"""
    create_filter_config(; kwargs...) -> FilterConfig

Create a FilterConfig instance with support for custom parameters

# Example
```julia
config = create_filter_config(
    min_total_depth = 30,
    max_eas_af = 0.01,
    min_revel_score = 0.75
)
```
"""
function create_filter_config(;
    min_total_depth::Int = 30,
    min_variant_frequency::Float64 = 0.03,
    max_eas_af::Float64 = 0.01,
    min_revel_score::Float64 = 0.75,
    min_primate_ai_score::Float64 = 0.8,
    min_dann_score::Float64 = 0.96
)
    config = FilterConfig(
        min_total_depth = min_total_depth,
        min_variant_frequency = min_variant_frequency,
        max_eas_af = max_eas_af,
        min_revel_score = min_revel_score,
        min_primate_ai_score = min_primate_ai_score,
        min_dann_score = min_dann_score
    )

    validate_config(config)
    return config
end

"""
    validate_config(config::FilterConfig)

Validate configuration parameters for reasonableness, throw error if unreasonable

# Validation rules
- All thresholds must be within reasonable ranges
- Frequencies must be between 0-1
- Depth must be positive integer
"""
function validate_config(config::FilterConfig)
    # Validate depth
    if config.min_total_depth < 1
        error("min_total_depth must be at least 1, got $(config.min_total_depth)")
    end

    # Validate VAF
    if config.min_variant_frequency < 0.0 || config.min_variant_frequency > 1.0
        error("min_variant_frequency must be between 0 and 1, got $(config.min_variant_frequency)")
    end

    # Validate population frequency
    if config.max_eas_af < 0.0 || config.max_eas_af > 1.0
        error("max_eas_af must be between 0 and 1, got $(config.max_eas_af)")
    end

    # Validate REVEL score
    if config.min_revel_score < 0.0 || config.min_revel_score > 1.0
        error("min_revel_score must be between 0 and 1, got $(config.min_revel_score)")
    end

    # Validate PrimateAI score
    if config.min_primate_ai_score < 0.0 || config.min_primate_ai_score > 1.0
        error("min_primate_ai_score must be between 0 and 1, got $(config.min_primate_ai_score)")
    end

    # Validate DANN score
    if config.min_dann_score < 0.0 || config.min_dann_score > 1.0
        error("min_dann_score must be between 0 and 1, got $(config.min_dann_score)")
    end

    return nothing
end

"""
    display_config(config::FilterConfig; io::IO = stdout)

Display configuration parameters in readable format
"""
function display_config(config::FilterConfig; io::IO = stdout)
    println(io, "=".^60)
    println(io, "JSON2MAF Filter Configuration")
    println(io, "=".^60)
    println(io)

    println(io, "Quality filtering parameters:")
    println(io, "  Minimum sequencing depth (min_total_depth):       $(config.min_total_depth)")
    println(io, "  Minimum VAF (min_variant_frequency):              $(config.min_variant_frequency)")
    println(io)

    println(io, "Population frequency filtering parameters:")
    println(io, "  Maximum East Asian AF (max_eas_af):               $(config.max_eas_af)")
    println(io)

    println(io, "Predictive score thresholds:")
    println(io, "  REVEL minimum score (min_revel_score):            $(config.min_revel_score)")
    println(io, "  PrimateAI-3D minimum score:                       $(config.min_primate_ai_score)")
    println(io, "  DANN minimum score:                               $(config.min_dann_score)")
    println(io)

    println(io, "=".^60)
end

end # module FilterConfigModule
