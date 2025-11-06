"""
ClinVarFilter.jl

ClinVar pathogenicity assessment and priority determination module
Determines variant pathogenicity based on review status and clinical significance
"""

module ClinVarFilter

using ..DataStructures
using Dates

export assess_clinvar_pathogenicity, get_review_status_priority,
       is_cancer_related, resolve_conflicting_entries

"""
    assess_clinvar_pathogenicity(entries::Vector{ClinVarEntry}) -> ClinVarAssessment

Assess pathogenicity of ClinVar annotations

# Assessment logic
1. Filter entries that are Pathogenic or Likely pathogenic
2. If multiple qualifying entries exist, select the highest priority entry
3. Priority based on: review status > cancer-related > last updated time

# Parameters
- `entries`: List of ClinVar annotations

# Returns
- `ClinVarAssessment`: Contains pathogenicity determination and selected entry

# Example
```julia
assessment = assess_clinvar_pathogenicity(variant.clinvar)
if assessment.is_pathogenic
    println("Variant is pathogenic")
end
```
"""
function assess_clinvar_pathogenicity(entries::Vector{ClinVarEntry})::ClinVarAssessment
    # If no ClinVar annotations, return negative result
    if isempty(entries)
        return ClinVarAssessment(
            false, false, nothing, "none",
            "No ClinVar entries available"
        )
    end

    # Filter Pathogenic and Likely pathogenic entries
    pathogenic_entries = filter(e -> is_pathogenic_entry(e), entries)

    # If no pathogenic entries
    if isempty(pathogenic_entries)
        return ClinVarAssessment(
            false, false, nothing, "none",
            "No pathogenic or likely pathogenic ClinVar entries"
        )
    end

    # Select best entry
    selected = resolve_conflicting_entries(pathogenic_entries)

    # Determine pathogenicity level
    # Strategy: If any form of "Pathogenic" is mentioned (including "Pathogenic/Likely pathogenic"),
    # consider as Pathogenic (adopting stronger evidence)
    # Only pure "Likely pathogenic" (without standalone "Pathogenic") is considered Likely pathogenic
    sig_lower = lowercase(selected.clinical_significance)

    # Check if contains standalone "Pathogenic" (not just "Likely pathogenic")
    # Split by delimiters and check each part
    parts = split(sig_lower, r"[/,;]")
    has_standalone_pathogenic = any(p -> occursin("pathogenic", strip(p)) && !occursin("likely", strip(p)), parts)
    has_likely_pathogenic = any(p -> occursin("likely", strip(p)) && occursin("pathogenic", strip(p)), parts)

    # If has standalone Pathogenic, classify as Pathogenic
    # Otherwise, if has Likely pathogenic, classify as Likely pathogenic
    is_path = has_standalone_pathogenic
    is_likely_path = !has_standalone_pathogenic && has_likely_pathogenic

    # Determine confidence level
    confidence = get_confidence_level(selected.review_status)

    # Build assessment reason
    reason = build_assessment_reason(selected, length(pathogenic_entries))

    return ClinVarAssessment(
        is_path,
        is_likely_path,
        selected,
        confidence,
        reason
    )
end

"""
    is_pathogenic_entry(entry::ClinVarEntry) -> Bool

Determine if ClinVar entry is pathogenic (Pathogenic or Likely pathogenic)
"""
function is_pathogenic_entry(entry::ClinVarEntry)::Bool
    if entry.clinical_significance === nothing
        return false
    end
    sig_lower = lowercase(entry.clinical_significance)
    return contains(sig_lower, "pathogenic") &&
           !contains(sig_lower, "benign") &&
           !contains(sig_lower, "uncertain")
end

"""
    get_review_status_priority(status::String) -> Int

Return review status priority score (lower number = higher priority)

# Review Status Priority Order
1. practice guideline
2. reviewed by expert panel
3. criteria provided, multiple submitters, no conflicts
4. criteria provided, conflicting interpretations
5. criteria provided, single submitter
6. no assertion criteria provided
7. no assertion provided
8. Other unknown status
"""
function get_review_status_priority(status::String)::Int
    status_lower = lowercase(status)

    if contains(status_lower, "practice guideline")
        return 1
    elseif contains(status_lower, "reviewed by expert panel")
        return 2
    elseif contains(status_lower, "multiple submitters") &&
           contains(status_lower, "no conflict")
        return 3
    elseif contains(status_lower, "conflicting")
        return 4
    elseif contains(status_lower, "single submitter")
        return 5
    elseif contains(status_lower, "no assertion criteria provided")
        return 6
    elseif contains(status_lower, "no assertion provided")
        return 7
    else
        return 8  # Unknown status, lowest priority
    end
end

"""
    is_cancer_related(diseases::Vector{String}) -> Bool

Determine if disease list contains cancer-related diseases

# Cancer keywords
- cancer, carcinoma, tumor, tumour, malignant, neoplasm
- lymphoma, leukemia, leukaemia, sarcoma, melanoma
- glioma, blastoma, myeloma
"""
function is_cancer_related(diseases::Vector{String})::Bool
    cancer_keywords = [
        "cancer", "carcinoma", "tumor", "tumour", "malignant", "neoplasm",
        "lymphoma", "leukemia", "leukaemia", "sarcoma", "melanoma",
        "glioma", "blastoma", "myeloma", "adenocarcinoma"
    ]

    for disease in diseases
        disease_lower = lowercase(disease)
        for keyword in cancer_keywords
            if contains(disease_lower, keyword)
                return true
            end
        end
    end

    return false
end

"""
    resolve_conflicting_entries(entries::Vector{ClinVarEntry}) -> ClinVarEntry

When multiple ClinVar annotations exist, select the best entry

# Priority rules
1. Review status priority highest (practice guideline > expert panel > ...)
2. Cancer-related diseases prioritized
3. Most recent update time prioritized
"""
function resolve_conflicting_entries(entries::Vector{ClinVarEntry})::ClinVarEntry
    if length(entries) == 1
        return entries[1]
    end

    # Sorting strategy:
    # 1. Review status priority (lower is better)
    # 2. Cancer-related (true > false)
    # 3. Last update time (newer is better)

    sorted = sort(entries, lt = (a, b) -> begin
        # First compare by review status priority (lower is better)
        pa = get_review_status_priority(a.review_status)
        pb = get_review_status_priority(b.review_status)
        if pa != pb
            return pa < pb
        end

        # Then compare by cancer-related (cancer-related comes first)
        ca = is_cancer_related(a.diseases)
        cb = is_cancer_related(b.diseases)
        if ca != cb
            return ca > cb
        end

        # Finally compare by last_evaluated (newer comes first)
        da = a.last_evaluated === nothing ? "" : a.last_evaluated
        db = b.last_evaluated === nothing ? "" : b.last_evaluated
        return da > db
    end)

    return sorted[1]
end

"""
    get_confidence_level(review_status::String) -> String

Determine confidence level based on review status
"""
function get_confidence_level(review_status::String)::String
    priority = get_review_status_priority(review_status)

    if priority <= 2
        return "high"
    elseif priority <= 4
        return "medium"
    else
        return "low"
    end
end

"""
    build_assessment_reason(entry::ClinVarEntry, total_entries::Int) -> String

Build assessment reason description
"""
function build_assessment_reason(entry::ClinVarEntry, total_entries::Int)::String
    parts = String[]

    # ClinVar classification
    push!(parts, "ClinVar: $(entry.clinical_significance)")

    # Review status
    push!(parts, "Review: $(entry.review_status)")

    # If multiple entries
    if total_entries > 1
        push!(parts, "Selected from $(total_entries) entries")
    end

    # Cancer-related
    if is_cancer_related(entry.diseases)
        push!(parts, "Cancer-related disease")
    end

    return join(parts, "; ")
end

end # module ClinVarFilter
