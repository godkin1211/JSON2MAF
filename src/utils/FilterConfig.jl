"""
FilterConfig.jl

過濾配置管理模組
提供配置的創建、驗證和顯示功能
"""

module FilterConfigModule

using ..DataStructures

export create_filter_config, validate_config, display_config

"""
    create_filter_config(; kwargs...) -> FilterConfig

創建一個 FilterConfig 實例，支援自訂參數

# 範例
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

驗證配置參數的合理性，如果不合理則拋出錯誤

# 驗證規則
- 所有閾值必須在合理範圍內
- 頻率必須在 0-1 之間
- 深度必須為正整數
"""
function validate_config(config::FilterConfig)
    # 驗證深度
    if config.min_total_depth < 1
        error("min_total_depth must be at least 1, got $(config.min_total_depth)")
    end

    # 驗證 VAF
    if config.min_variant_frequency < 0.0 || config.min_variant_frequency > 1.0
        error("min_variant_frequency must be between 0 and 1, got $(config.min_variant_frequency)")
    end

    # 驗證族群頻率
    if config.max_eas_af < 0.0 || config.max_eas_af > 1.0
        error("max_eas_af must be between 0 and 1, got $(config.max_eas_af)")
    end

    # 驗證 REVEL 分數
    if config.min_revel_score < 0.0 || config.min_revel_score > 1.0
        error("min_revel_score must be between 0 and 1, got $(config.min_revel_score)")
    end

    # 驗證 PrimateAI 分數
    if config.min_primate_ai_score < 0.0 || config.min_primate_ai_score > 1.0
        error("min_primate_ai_score must be between 0 and 1, got $(config.min_primate_ai_score)")
    end

    # 驗證 DANN 分數
    if config.min_dann_score < 0.0 || config.min_dann_score > 1.0
        error("min_dann_score must be between 0 and 1, got $(config.min_dann_score)")
    end

    return nothing
end

"""
    display_config(config::FilterConfig; io::IO = stdout)

以易讀格式顯示配置參數
"""
function display_config(config::FilterConfig; io::IO = stdout)
    println(io, "=".^60)
    println(io, "JSON2MAF Filter Configuration")
    println(io, "=".^60)
    println(io)

    println(io, "品質過濾參數:")
    println(io, "  最小測序深度 (min_total_depth):       $(config.min_total_depth)")
    println(io, "  最小VAF (min_variant_frequency):      $(config.min_variant_frequency)")
    println(io)

    println(io, "族群頻率過濾參數:")
    println(io, "  東亞族群最大AF (max_eas_af):          $(config.max_eas_af)")
    println(io)

    println(io, "預測分數閾值:")
    println(io, "  REVEL 最小分數 (min_revel_score):     $(config.min_revel_score)")
    println(io, "  PrimateAI-3D 最小分數:                $(config.min_primate_ai_score)")
    println(io, "  DANN 最小分數:                        $(config.min_dann_score)")
    println(io)

    println(io, "=".^60)
end

end # module FilterConfigModule
