"""
DecisionEngine.jl

Integrated filtering decision engine
Integrates ClinVar and predictive score assessment to make final inclusion/exclusion decisions
"""

module DecisionEngine

using ..DataStructures

export make_filter_decision, count_supporting_predictive_scores, has_primate_ai_support

"""
    make_filter_decision(variant::VariantPosition,
                        clinvar_assessment::ClinVarAssessment,
                        predictive_assessment::PredictiveAssessment) -> FilterDecision

Integrate ClinVar and predictive score assessment to make final filtering decision

# Decision logic
1. ClinVar Pathogenic → Include directly, mark as "Pathogenic"
2. ClinVar Likely pathogenic → Include directly, mark as "Likely pathogenic"
3. ClinVar inconclusive + 2+ predictive scores support → Include, mark as "Likely pathogenic"
4. ClinVar inconclusive + only PrimateAI-3D support → Include, mark as "Likely pathogenic"
5. Other cases → Exclude

# Parameters
- `variant`: Variant position data
- `clinvar_assessment`: ClinVar assessment result
- `predictive_assessment`: Predictive score assessment result

# Returns
- `FilterDecision`: Contains whether to include, pathogenicity classification, evidence source, justification

# Example
```julia
decision = make_filter_decision(variant, clinvar_result, predictive_result)
if decision.should_include
    println("Include variant: ", decision.pathogenicity_class)
    println("Reason: ", decision.justification)
end
```
"""
function make_filter_decision(
    variant::VariantPosition,
    clinvar_assessment::ClinVarAssessment,
    predictive_assessment::PredictiveAssessment
)::FilterDecision

    # Rule 1: ClinVar Pathogenic
    if clinvar_assessment.is_pathogenic
        return FilterDecision(
            true,
            "Pathogenic",
            "ClinVar",
            build_clinvar_justification(clinvar_assessment, "Pathogenic")
        )
    end

    # Rule 2: ClinVar Likely pathogenic
    if clinvar_assessment.is_likely_pathogenic
        return FilterDecision(
            true,
            "Likely pathogenic",
            "ClinVar",
            build_clinvar_justification(clinvar_assessment, "Likely pathogenic")
        )
    end

    # ClinVar inconclusive, check predictive scores
    # Rules 3 & 4: Predictive score support
    if predictive_assessment.suggests_pathogenic
        # Check if PrimateAI-3D alone supports
        if has_primate_ai_support(predictive_assessment) &&
           predictive_assessment.support_count == 1
            return FilterDecision(
                true,
                "Likely pathogenic",
                "Predictive",
                "Supported by PrimateAI-3D (highest priority predictor)"
            )
        end

        # 2+ predictive scores support
        if predictive_assessment.support_count >= 2
            return FilterDecision(
                true,
                "Likely pathogenic",
                "Predictive",
                build_predictive_justification(predictive_assessment)
            )
        end
    end

    # Rule 5: Other cases → Exclude
    return FilterDecision(
        false,
        "Excluded",
        "Insufficient evidence",
        build_exclusion_justification(clinvar_assessment, predictive_assessment)
    )
end

"""
    count_supporting_predictive_scores(assessment::PredictiveAssessment) -> Int

Count number of predictive scores supporting pathogenicity
"""
function count_supporting_predictive_scores(assessment::PredictiveAssessment)::Int
    return assessment.support_count
end

"""
    has_primate_ai_support(assessment::PredictiveAssessment) -> Bool

Determine if PrimateAI-3D support exists
"""
function has_primate_ai_support(assessment::PredictiveAssessment)::Bool
    return assessment.has_primate_ai_support
end

"""
    build_clinvar_justification(assessment::ClinVarAssessment, classification::String) -> String

Build justification description for ClinVar evidence
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

Build justification description for predictive score evidence
"""
function build_predictive_justification(assessment::PredictiveAssessment)::String
    parts = String[]

    push!(parts, "Supported by $(assessment.support_count) predictive scores")

    # List supporting scores
    score_names = sort(collect(keys(assessment.contributing_scores)))
    push!(parts, "Scores: $(join(score_names, ", "))")

    push!(parts, "Confidence: $(round(assessment.confidence, digits=2))")

    return join(parts, "; ")
end

"""
    build_exclusion_justification(clinvar_assessment::ClinVarAssessment,
                                  predictive_assessment::PredictiveAssessment) -> String

Build justification description for excluded variant
"""
function build_exclusion_justification(
    clinvar_assessment::ClinVarAssessment,
    predictive_assessment::PredictiveAssessment
)::String
    reasons = String[]

    # ClinVar status
    if !clinvar_assessment.is_pathogenic && !clinvar_assessment.is_likely_pathogenic
        push!(reasons, "No pathogenic ClinVar evidence")
    end

    # Predictive score status
    support_count = predictive_assessment.support_count
    if support_count == 0
        push!(reasons, "No supporting predictive scores")
    elseif support_count == 1
        # Check which score
        if has_primate_ai_support(predictive_assessment)
            # Should not reach here, as PrimateAI-3D alone should pass
            push!(reasons, "Only PrimateAI-3D support (should not reach here)")
        else
            score_name = first(keys(predictive_assessment.contributing_scores))
            push!(reasons, "Only $(score_name) support (insufficient, need 2+ or PrimateAI-3D)")
        end
    end

    return join(reasons, "; ")
end

end # module DecisionEngine
