"""
DecisionEngine.jl

整合過濾決策引擎
整合 ClinVar 與預測分數評估，做出最終納入/排除決策
"""

module DecisionEngine

using ..DataStructures

export make_filter_decision, count_supporting_predictive_scores, has_primate_ai_support

"""
    make_filter_decision(variant::VariantPosition,
                        clinvar_assessment::ClinVarAssessment,
                        predictive_assessment::PredictiveAssessment) -> FilterDecision

整合 ClinVar 與預測分數評估，做出最終過濾決策

# 決策邏輯
1. ClinVar Pathogenic → 直接納入，標記為 "Pathogenic"
2. ClinVar Likely pathogenic → 直接納入，標記為 "Likely pathogenic"
3. ClinVar 無結論 + 2+ 預測分數支持 → 納入，標記為 "Likely pathogenic"
4. ClinVar 無結論 + 僅 PrimateAI-3D 支持 → 納入，標記為 "Likely pathogenic"
5. 其他情況 → 排除

# 參數
- `variant`: 變異位點資料
- `clinvar_assessment`: ClinVar 評估結果
- `predictive_assessment`: 預測分數評估結果

# 返回
- `FilterDecision`: 包含是否納入、致病性分類、證據來源、理由

# 範例
```julia
decision = make_filter_decision(variant, clinvar_result, predictive_result)
if decision.should_include
    println("納入變異: ", decision.pathogenicity_class)
    println("理由: ", decision.justification)
end
```
"""
function make_filter_decision(
    variant::VariantPosition,
    clinvar_assessment::ClinVarAssessment,
    predictive_assessment::PredictiveAssessment
)::FilterDecision

    # 規則 1: ClinVar Pathogenic
    if clinvar_assessment.is_pathogenic
        return FilterDecision(
            true,
            "Pathogenic",
            "ClinVar",
            build_clinvar_justification(clinvar_assessment, "Pathogenic")
        )
    end

    # 規則 2: ClinVar Likely pathogenic
    if clinvar_assessment.is_likely_pathogenic
        return FilterDecision(
            true,
            "Likely pathogenic",
            "ClinVar",
            build_clinvar_justification(clinvar_assessment, "Likely pathogenic")
        )
    end

    # ClinVar 無結論，檢查預測分數
    # 規則 3 & 4: 預測分數支持
    if predictive_assessment.suggests_pathogenic
        # 檢查是否為 PrimateAI-3D 單獨支持
        if has_primate_ai_support(predictive_assessment) &&
           predictive_assessment.support_count == 1
            return FilterDecision(
                true,
                "Likely pathogenic",
                "Predictive",
                "Supported by PrimateAI-3D (highest priority predictor)"
            )
        end

        # 2+ 預測分數支持
        if predictive_assessment.support_count >= 2
            return FilterDecision(
                true,
                "Likely pathogenic",
                "Predictive",
                build_predictive_justification(predictive_assessment)
            )
        end
    end

    # 規則 5: 其他情況 → 排除
    return FilterDecision(
        false,
        "Excluded",
        "Insufficient evidence",
        build_exclusion_justification(clinvar_assessment, predictive_assessment)
    )
end

"""
    count_supporting_predictive_scores(assessment::PredictiveAssessment) -> Int

計算支持致病性的預測分數數量
"""
function count_supporting_predictive_scores(assessment::PredictiveAssessment)::Int
    return assessment.support_count
end

"""
    has_primate_ai_support(assessment::PredictiveAssessment) -> Bool

判斷是否有 PrimateAI-3D 支持
"""
function has_primate_ai_support(assessment::PredictiveAssessment)::Bool
    return assessment.has_primate_ai_support
end

"""
    build_clinvar_justification(assessment::ClinVarAssessment, classification::String) -> String

建立 ClinVar 證據的理由說明
"""
function build_clinvar_justification(assessment::ClinVarAssessment, classification::String)::String
    parts = String[]

    push!(parts, "ClinVar: $(classification)")

    if assessment.selected_entry !== nothing
        push!(parts, "Review status: $(assessment.selected_entry.review_status)")
        push!(parts, "Confidence: $(assessment.confidence_level)")
    end

    return join(parts, "; ")
end

"""
    build_predictive_justification(assessment::PredictiveAssessment) -> String

建立預測分數證據的理由說明
"""
function build_predictive_justification(assessment::PredictiveAssessment)::String
    parts = String[]

    push!(parts, "Supported by $(assessment.support_count) predictive scores")

    # 列出支持的分數
    score_names = sort(collect(keys(assessment.contributing_scores)))
    push!(parts, "Scores: $(join(score_names, ", "))")

    push!(parts, "Confidence: $(round(assessment.confidence, digits=2))")

    return join(parts, "; ")
end

"""
    build_exclusion_justification(clinvar_assessment::ClinVarAssessment,
                                  predictive_assessment::PredictiveAssessment) -> String

建立排除變異的理由說明
"""
function build_exclusion_justification(
    clinvar_assessment::ClinVarAssessment,
    predictive_assessment::PredictiveAssessment
)::String
    reasons = String[]

    # ClinVar 狀態
    if !clinvar_assessment.is_pathogenic && !clinvar_assessment.is_likely_pathogenic
        push!(reasons, "No pathogenic ClinVar evidence")
    end

    # 預測分數狀態
    support_count = predictive_assessment.support_count
    if support_count == 0
        push!(reasons, "No supporting predictive scores")
    elseif support_count == 1
        # 檢查是哪個分數
        if has_primate_ai_support(predictive_assessment)
            # 不應該到這裡，因為 PrimateAI-3D 單獨就應該通過
            push!(reasons, "Only PrimateAI-3D support (should not reach here)")
        else
            score_name = first(keys(predictive_assessment.contributing_scores))
            push!(reasons, "Only $(score_name) support (insufficient, need 2+ or PrimateAI-3D)")
        end
    end

    return join(reasons, "; ")
end

end # module DecisionEngine
