"""
PredictiveScores.jl

Predictive score assessment module
Integrates multiple predictive scores (PrimateAI-3D, DANN, REVEL, COSMIC)
Provides supplementary determination for likely pathogenic
"""

module PredictiveScores

using ..DataStructures

export assess_predictive_scores, get_primate_ai_score, get_dann_score,
       get_revel_score, is_in_cosmic

"""
    assess_predictive_scores(variant::VariantPosition, config::FilterConfig) -> PredictiveAssessment

Assess if predictive scores support pathogenicity

# Assessment logic
1. Check if each predictive score exceeds threshold
2. Count number of scores supporting pathogenicity
3. Determine if should be suggested as likely pathogenic

# Support conditions (any one met suggests likely pathogenic)
- PrimateAI-3D ≥ threshold (highest priority, alone sufficient)
- 2+ other scores ≥ threshold (REVEL, DANN, COSMIC)

# Parameters
- `variant`: Variant position data
- `config`: Filter configuration parameters (including predictive score thresholds)

# Returns
- `PredictiveAssessment`: Contains whether pathogenicity is suggested and contributing scores

# Example
```julia
assessment = assess_predictive_scores(variant, config)
if assessment.suggests_pathogenic
    println("Predictive scores support pathogenicity")
    println("Contributing scores: ", assessment.contributing_scores)
end
```
"""
function assess_predictive_scores(variant::VariantPosition, config::FilterConfig)::PredictiveAssessment
    contributing = Dict{String, Float64}()
    support_count = 0

    # Check PrimateAI-3D
    primate_ai_3d = get_primate_ai_score(variant)
    if primate_ai_3d !== nothing && primate_ai_3d >= config.min_primate_ai_score
        contributing["PrimateAI-3D"] = primate_ai_3d
        support_count += 1
    end

    # Check REVEL
    revel = get_revel_score(variant)
    if revel !== nothing && revel >= config.min_revel_score
        contributing["REVEL"] = revel
        support_count += 1
    end

    # Check DANN
    dann = get_dann_score(variant)
    if dann !== nothing && dann >= config.min_dann_score
        contributing["DANN"] = dann
        support_count += 1
    end

    # Check COSMIC (presence indicates positive evidence)
    if is_in_cosmic(variant)
        contributing["COSMIC"] = 1.0  # Use 1.0 to indicate presence
        support_count += 1
    end

    # Determine if should be suggested as likely pathogenic
    # Condition 1: PrimateAI-3D alone supports (highest priority)
    # Condition 2: 2+ other scores support
    has_primate_ai_3d = haskey(contributing, "PrimateAI-3D")
    suggests_pathogenic = has_primate_ai_3d || support_count >= 2

    # Calculate confidence score (0-1)
    # Based on number and types of supporting scores
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

Calculate confidence score for predictive assessment

# Confidence levels
- PrimateAI-3D present: base 0.7
- Each additional supporting score: +0.1
- Maximum: 1.0
"""
function calculate_confidence(contributing::Dict{String, Float64}, has_primate_ai::Bool, support_count::Int)::Float64
    if support_count == 0
        return 0.0
    end

    # If PrimateAI-3D present, base confidence is 0.7
    base = has_primate_ai ? 0.7 : 0.5

    # Each additional supporting score increases confidence
    confidence = base + (support_count - 1) * 0.1

    # Maximum 1.0
    return min(confidence, 1.0)
end

"""
    get_primate_ai_score(variant::VariantPosition) -> Union{Float64, Nothing}

Extract PrimateAI-3D score (prioritize 3D version)
"""
function get_primate_ai_score(variant::VariantPosition)::Union{Float64, Nothing}
    # Prioritize PrimateAI-3D
    if variant.primate_ai_3d !== nothing
        return variant.primate_ai_3d
    end

    # Fallback: use PrimateAI (legacy version)
    if variant.primate_ai !== nothing
        return variant.primate_ai
    end

    return nothing
end

"""
    get_dann_score(variant::VariantPosition) -> Union{Float64, Nothing}

Extract DANN score
"""
function get_dann_score(variant::VariantPosition)::Union{Float64, Nothing}
    return variant.dann_score
end

"""
    get_revel_score(variant::VariantPosition) -> Union{Float64, Nothing}

Extract REVEL score
"""
function get_revel_score(variant::VariantPosition)::Union{Float64, Nothing}
    return variant.revel_score
end

"""
    is_in_cosmic(variant::VariantPosition) -> Bool

Determine if variant is in COSMIC database
"""
function is_in_cosmic(variant::VariantPosition)::Bool
    return !isempty(variant.cosmic)
end

end # module PredictiveScores
