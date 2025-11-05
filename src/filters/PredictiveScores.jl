"""
PredictiveScores.jl

預測分數評估模組
整合多種預測分數（PrimateAI-3D, DANN, REVEL, COSMIC）
提供 likely pathogenic 補充判斷
"""

module PredictiveScores

using ..DataStructures

export assess_predictive_scores, get_primate_ai_score, get_dann_score,
       get_revel_score, is_in_cosmic

"""
    assess_predictive_scores(variant::VariantPosition, config::FilterConfig) -> PredictiveAssessment

評估預測分數是否支持致病性

# 評估邏輯
1. 檢查各項預測分數是否超過閾值
2. 統計支持致病性的分數數量
3. 判斷是否建議為 likely pathogenic

# 支持條件（任一符合即建議為 likely pathogenic）
- PrimateAI-3D ≥ 閾值（最高優先級，單獨即可）
- 2+ 其他分數 ≥ 閾值（REVEL, DANN, COSMIC）

# 參數
- `variant`: 變異位點資料
- `config`: 過濾配置參數（包含預測分數閾值）

# 返回
- `PredictiveAssessment`: 包含是否建議致病性及貢獻分數

# 範例
```julia
assessment = assess_predictive_scores(variant, config)
if assessment.suggests_pathogenic
    println("預測分數支持致病性")
    println("貢獻分數: ", assessment.contributing_scores)
end
```
"""
function assess_predictive_scores(variant::VariantPosition, config::FilterConfig)::PredictiveAssessment
    contributing = Dict{String, Float64}()
    support_count = 0

    # 檢查 PrimateAI-3D
    primate_ai_3d = get_primate_ai_score(variant)
    if primate_ai_3d !== nothing && primate_ai_3d >= config.min_primate_ai_score
        contributing["PrimateAI-3D"] = primate_ai_3d
        support_count += 1
    end

    # 檢查 REVEL
    revel = get_revel_score(variant)
    if revel !== nothing && revel >= config.min_revel_score
        contributing["REVEL"] = revel
        support_count += 1
    end

    # 檢查 DANN
    dann = get_dann_score(variant)
    if dann !== nothing && dann >= config.min_dann_score
        contributing["DANN"] = dann
        support_count += 1
    end

    # 檢查 COSMIC（存在即為陽性證據）
    if is_in_cosmic(variant)
        contributing["COSMIC"] = 1.0  # 用 1.0 表示存在
        support_count += 1
    end

    # 判斷是否建議為 likely pathogenic
    # 條件1：PrimateAI-3D 單獨支持（最高優先級）
    # 條件2：2+ 其他分數支持
    has_primate_ai_3d = haskey(contributing, "PrimateAI-3D")
    suggests_pathogenic = has_primate_ai_3d || support_count >= 2

    # 計算信心分數（0-1）
    # 基於支持的分數數量和種類
    confidence = calculate_confidence(contributing, has_primate_ai_3d, support_count)

    return PredictiveAssessment(
        suggests_pathogenic,
        contributing,
        confidence,
        support_count,
        has_primate_ai_3d
    )
end

"""
    calculate_confidence(contributing::Dict{String, Float64}, has_primate_ai::Bool, support_count::Int) -> Float64

計算預測評估的信心分數

# 信心等級
- PrimateAI-3D 存在: 基礎 0.7
- 每多一個支持分數: +0.1
- 最高: 1.0
"""
function calculate_confidence(contributing::Dict{String, Float64}, has_primate_ai::Bool, support_count::Int)::Float64
    if support_count == 0
        return 0.0
    end

    # 如果有 PrimateAI-3D，基礎信心為 0.7
    base = has_primate_ai ? 0.7 : 0.5

    # 每增加一個支持分數，信心增加
    confidence = base + (support_count - 1) * 0.1

    # 最高 1.0
    return min(confidence, 1.0)
end

"""
    get_primate_ai_score(variant::VariantPosition) -> Union{Float64, Nothing}

提取 PrimateAI-3D 分數（優先使用 3D 版本）
"""
function get_primate_ai_score(variant::VariantPosition)::Union{Float64, Nothing}
    # 優先使用 PrimateAI-3D
    if variant.primate_ai_3d !== nothing
        return variant.primate_ai_3d
    end

    # 備用：使用 PrimateAI（舊版）
    if variant.primate_ai !== nothing
        return variant.primate_ai
    end

    return nothing
end

"""
    get_dann_score(variant::VariantPosition) -> Union{Float64, Nothing}

提取 DANN 分數
"""
function get_dann_score(variant::VariantPosition)::Union{Float64, Nothing}
    return variant.dann_score
end

"""
    get_revel_score(variant::VariantPosition) -> Union{Float64, Nothing}

提取 REVEL 分數
"""
function get_revel_score(variant::VariantPosition)::Union{Float64, Nothing}
    return variant.revel_score
end

"""
    is_in_cosmic(variant::VariantPosition) -> Bool

判斷變異是否在 COSMIC 資料庫中
"""
function is_in_cosmic(variant::VariantPosition)::Bool
    return !isempty(variant.cosmic)
end

"""
    format_contributing_scores(contributing::Dict{String, Float64}) -> String

格式化貢獻分數為字串（用於報告）
"""
function format_contributing_scores(contributing::Dict{String, Float64})::String
    if isempty(contributing)
        return "None"
    end

    parts = String[]
    for (name, score) in sort(collect(contributing), by = x -> x[1])
        if name == "COSMIC"
            push!(parts, "COSMIC: present")
        else
            push!(parts, "$(name): $(round(score, digits=3))")
        end
    end

    return join(parts, ", ")
end

end # module PredictiveScores
