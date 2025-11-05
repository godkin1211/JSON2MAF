"""
ClinVarFilter.jl

ClinVar 致病性評估與優先級判斷模組
根據 review status 與 clinical significance 判斷變異的致病性
"""

module ClinVarFilter

using ..DataStructures
using Dates

export assess_clinvar_pathogenicity, get_review_status_priority,
       is_cancer_related, resolve_conflicting_entries

"""
    assess_clinvar_pathogenicity(entries::Vector{ClinVarEntry}) -> ClinVarAssessment

評估 ClinVar 註解的致病性

# 評估邏輯
1. 過濾出 Pathogenic 或 Likely pathogenic 的條目
2. 如有多筆符合條目，選擇最高優先級的條目
3. 優先級判斷依據：review status > 癌症相關 > 最後更新時間

# 參數
- `entries`: ClinVar 註解列表

# 返回
- `ClinVarAssessment`: 包含致病性判斷與選定條目

# 範例
```julia
assessment = assess_clinvar_pathogenicity(variant.clinvar)
if assessment.is_pathogenic
    println("變異為致病性")
end
```
"""
function assess_clinvar_pathogenicity(entries::Vector{ClinVarEntry})::ClinVarAssessment
    # 如果沒有 ClinVar 註解，返回陰性結果
    if isempty(entries)
        return ClinVarAssessment(
            false, false, nothing, "none",
            "No ClinVar entries available"
        )
    end

    # 過濾 Pathogenic 和 Likely pathogenic 條目
    pathogenic_entries = filter(e -> is_pathogenic_entry(e), entries)

    # 如果沒有致病性條目
    if isempty(pathogenic_entries)
        return ClinVarAssessment(
            false, false, nothing, "none",
            "No pathogenic or likely pathogenic ClinVar entries"
        )
    end

    # 選擇最佳條目
    selected = resolve_conflicting_entries(pathogenic_entries)

    # 判斷致病性等級
    # 策略：如果提到任何形式的 "Pathogenic"（包括 "Pathogenic/Likely pathogenic"），
    # 視為 Pathogenic（採用較強的證據）
    # 只有純粹的 "Likely pathogenic"（沒有單獨的 "Pathogenic"）才視為 Likely pathogenic
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

    # 判斷信心等級
    confidence = get_confidence_level(selected.review_status)

    # 建立理由說明
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

判斷 ClinVar 條目是否為致病性 (Pathogenic 或 Likely pathogenic)
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

返回 review status 的優先級分數 (數字越小優先級越高)

# Review Status 優先級排序
1. practice guideline
2. reviewed by expert panel
3. criteria provided, multiple submitters, no conflicts
4. criteria provided, conflicting interpretations
5. criteria provided, single submitter
6. no assertion criteria provided
7. no assertion provided
8. 其他未知狀態
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
        return 8  # 未知狀態，最低優先級
    end
end

"""
    is_cancer_related(diseases::Vector{String}) -> Bool

判斷疾病列表中是否包含癌症相關疾病

# 癌症關鍵字
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

當多筆 ClinVar 註解存在時，選擇最佳條目

# 優先級規則
1. Review status 優先級最高 (practice guideline > expert panel > ...)
2. 癌症相關疾病優先
3. 最近更新時間優先
"""
function resolve_conflicting_entries(entries::Vector{ClinVarEntry})::ClinVarEntry
    if length(entries) == 1
        return entries[1]
    end

    # 排序策略：
    # 1. Review status 優先級 (越小越好)
    # 2. 癌症相關 (true > false)
    # 3. 最後更新時間 (越新越好)

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

根據 review status 判斷信心等級
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

建立評估理由說明
"""
function build_assessment_reason(entry::ClinVarEntry, total_entries::Int)::String
    parts = String[]

    # ClinVar 分類
    push!(parts, "ClinVar: $(entry.clinical_significance)")

    # Review status
    push!(parts, "Review: $(entry.review_status)")

    # 如有多筆條目
    if total_entries > 1
        push!(parts, "Selected from $(total_entries) entries")
    end

    # 癌症相關
    if is_cancer_related(entry.diseases)
        push!(parts, "Cancer-related disease")
    end

    return join(parts, "; ")
end

end # module ClinVarFilter
