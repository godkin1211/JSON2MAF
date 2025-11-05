"""
MAFConverter.jl

MAF 格式轉換模組
將 Nirvana 註解的變異轉換為 MAF (Mutation Annotation Format) 格式
"""

module MAFConverter

using ..DataStructures
import Base: get

export variant_to_maf, map_variant_classification, map_variant_type,
       extract_hgvs_notation, select_canonical_transcript

# 全局常數：避免重複分配空字串
const EMPTY_STRING = ""
const DOT_STRING = "."

# Consequence 嚴重程度映射（全局常數，避免每次函數調用都重新創建 Dict）
# 數字越小越嚴重，基於 Sequence Ontology 和 VEP 的分類
const CONSEQUENCE_SEVERITY = Dict{String, Int}(
    # HIGH impact
    "transcript_ablation" => 1,
    "splice_acceptor_variant" => 2,
    "splice_donor_variant" => 3,
    "stop_gained" => 4,
    "frameshift_variant" => 5,
    "stop_lost" => 6,
    "start_lost" => 7,
    "transcript_amplification" => 8,
    # MODERATE impact
    "inframe_insertion" => 10,
    "inframe_deletion" => 11,
    "missense_variant" => 12,
    "protein_altering_variant" => 13,
    # LOW impact
    "splice_region_variant" => 20,
    "incomplete_terminal_codon_variant" => 21,
    "start_retained_variant" => 22,
    "stop_retained_variant" => 23,
    "synonymous_variant" => 24,
    # MODIFIER
    "coding_sequence_variant" => 30,
    "mature_miRNA_variant" => 31,
    "5_prime_UTR_variant" => 32,
    "3_prime_UTR_variant" => 33,
    "non_coding_transcript_exon_variant" => 34,
    "intron_variant" => 35,
    "NMD_transcript_variant" => 36,
    "non_coding_transcript_variant" => 37,
    "upstream_gene_variant" => 40,
    "downstream_gene_variant" => 41,
    "TFBS_ablation" => 42,
    "TFBS_amplification" => 43,
    "TF_binding_site_variant" => 44,
    "regulatory_region_ablation" => 45,
    "regulatory_region_amplification" => 46,
    "feature_elongation" => 47,
    "regulatory_region_variant" => 48,
    "feature_truncation" => 49,
    "intergenic_variant" => 50
)

"""
    variant_to_maf(variant::VariantPosition, decision::FilterDecision) -> MAFRecord

將變異位點轉換為 MAF 記錄

# 參數
- `variant`: 變異位點資料
- `decision`: 過濾決策結果

# 返回
- `MAFRecord`: MAF 格式記錄

# 範例
```julia
maf_record = variant_to_maf(variant, decision)
```
"""
function variant_to_maf(variant::VariantPosition, decision::FilterDecision)::MAFRecord
    # 選擇 canonical transcript
    canonical = select_canonical_transcript(variant.transcripts)

    # 提取基因符號
    hugo_symbol = canonical !== nothing ? get(canonical, :gene_symbol, ".") : "."

    # 提取 HGVS
    hgvsc, hgvsp, hgvsp_short = if canonical !== nothing
        extract_hgvs_notation(canonical)
    else
        (".", ".", ".")
    end

    # 提取 Transcript ID
    transcript_id = canonical !== nothing ? get(canonical, :transcript_id, ".") : "."

    # 映射變異分類
    variant_classification = if canonical !== nothing
        consequences = get(canonical, :consequences, String[])
        if !isempty(consequences)
            map_variant_classification(consequences[1])
        else
            "."
        end
    else
        "."
    end

    # 映射變異類型
    variant_type_maf = map_variant_type(variant.reference_allele, variant.alternate_allele)

    # 提取 ClinVar 資訊
    clinvar_id, clinvar_review, clinvar_sig, clinvar_disease = extract_clinvar_info(variant.clinvar)

    # 提取 COSMIC ID
    cosmic_id = extract_cosmic_id(variant.cosmic)

    # 提取 dbSNP
    dbsnp_rs = isempty(variant.dbsnp_ids) ? "." : variant.dbsnp_ids[1]

    # 提取 gnomAD AF
    gnomad_af = extract_gnomad_af(variant.population_frequencies)

    # 提取預測分數並轉為字串 (使用 nothing 表示缺失值)
    # PrimateAI: 優先使用 3D 版本
    primate_ai = variant.primate_ai_3d !== nothing ? variant.primate_ai_3d :
                 (variant.primate_ai !== nothing ? variant.primate_ai : nothing)
    primate_ai_str = isnothing(primate_ai) ? nothing : string(primate_ai)
    dann_str = isnothing(variant.dann_score) ? nothing : string(variant.dann_score)
    revel_str = isnothing(variant.revel_score) ? nothing : string(variant.revel_score)
    gnomad_af_str = isnothing(gnomad_af) ? nothing : string(gnomad_af)

    # 提取東亞 AF
    gnomad_eas_af = extract_gnomad_eas_af(variant.population_frequencies)
    gnomad_eas_af_str = isnothing(gnomad_eas_af) ? nothing : string(gnomad_eas_af)

    # 提取深度與 VAF - 使用全局常數減少分配
    depth_str = isnothing(variant.total_depth) ? EMPTY_STRING : string(variant.total_depth)
    vaf_str = (isnothing(variant.variant_frequencies) || isempty(variant.variant_frequencies)) ?
              EMPTY_STRING : string(variant.variant_frequencies[1])

    return MAFRecord(
        hugo_symbol = hugo_symbol,
        chromosome = variant.chromosome,
        start_position = variant.start,
        end_position = variant.end_pos,
        strand = "+",  # Nirvana 通常不提供，預設為 +
        variant_classification = variant_classification,
        variant_type = variant_type_maf,
        reference_allele = variant.reference_allele,
        tumor_seq_allele1 = variant.reference_allele,  # 通常與 reference 相同
        tumor_seq_allele2 = variant.alternate_allele,  # 變異等位基因
        tumor_sample_barcode = "",  # 需從 samples 提取，目前設為空
        hgvsc = hgvsc,
        hgvsp = hgvsp,
        hgvsp_short = hgvsp_short,
        transcript_id = transcript_id,
        dbsnp_rs = dbsnp_rs,
        dbsnp_val_status = "",  # 不可用
        cosmic_id = cosmic_id,
        clinvar_id = clinvar_id,
        clinvar_review_status = clinvar_review,
        clinvar_significance = clinvar_sig,
        clinvar_disease = clinvar_disease,
        primate_ai_score = primate_ai_str,
        dann_score = dann_str,
        revel_score = revel_str,
        gnomad_af = gnomad_af_str,
        gnomad_eas_af = gnomad_eas_af_str,
        depth = depth_str,
        vaf = vaf_str
    )
end

"""
    map_variant_classification(consequence::String) -> String

映射 Nirvana consequence 到 MAF Variant_Classification

# MAF Variant Classification
- Missense_Mutation
- Nonsense_Mutation
- Frame_Shift_Del
- Frame_Shift_Ins
- In_Frame_Del
- In_Frame_Ins
- Splice_Site
- Translation_Start_Site
- Nonstop_Mutation
- Silent
- Intron
- RNA
- 3'UTR
- 5'UTR
- IGR (Intergenic_Region)
- etc.
"""
function map_variant_classification(consequence::String)::String
    consequence_lower = lowercase(consequence)

    # Missense
    if contains(consequence_lower, "missense")
        return "Missense_Mutation"

    # Nonsense (stop gained)
    elseif contains(consequence_lower, "stop_gained") || contains(consequence_lower, "nonsense")
        return "Nonsense_Mutation"

    # Frameshift
    elseif contains(consequence_lower, "frameshift")
        if contains(consequence_lower, "deletion") || contains(consequence_lower, "del")
            return "Frame_Shift_Del"
        elseif contains(consequence_lower, "insertion") || contains(consequence_lower, "ins")
            return "Frame_Shift_Ins"
        else
            return "Frame_Shift_Indel"
        end

    # Inframe indels
    elseif contains(consequence_lower, "inframe") || contains(consequence_lower, "in_frame")
        if contains(consequence_lower, "deletion") || contains(consequence_lower, "del")
            return "In_Frame_Del"
        elseif contains(consequence_lower, "insertion") || contains(consequence_lower, "ins")
            return "In_Frame_Ins"
        else
            return "In_Frame_Indel"
        end

    # Splice site
    elseif contains(consequence_lower, "splice")
        return "Splice_Site"

    # Start codon
    elseif contains(consequence_lower, "start") && contains(consequence_lower, "lost")
        return "Translation_Start_Site"

    # Stop lost
    elseif contains(consequence_lower, "stop") && contains(consequence_lower, "lost")
        return "Nonstop_Mutation"

    # Synonymous
    elseif contains(consequence_lower, "synonymous")
        return "Silent"

    # UTR
    elseif contains(consequence_lower, "3") && contains(consequence_lower, "utr")
        return "3'UTR"
    elseif contains(consequence_lower, "5") && contains(consequence_lower, "utr")
        return "5'UTR"

    # Intron
    elseif contains(consequence_lower, "intron")
        return "Intron"

    # Intergenic
    elseif contains(consequence_lower, "intergenic")
        return "IGR"

    # RNA
    elseif contains(consequence_lower, "rna") || contains(consequence_lower, "non_coding")
        return "RNA"

    # Default
    else
        return consequence  # 返回原始值
    end
end

"""
    map_variant_type(ref::String, alt::String) -> String

映射變異類型到 MAF Variant_Type

# MAF Variant Types
- SNP (single nucleotide polymorphism)
- DNP (di-nucleotide polymorphism)
- TNP (tri-nucleotide polymorphism)
- ONP (oligo-nucleotide polymorphism)
- INS (insertion)
- DEL (deletion)
- Consolidated (complex)
"""
function map_variant_type(ref::String, alt::String)::String
    ref_len = length(ref)
    alt_len = length(alt)

    # Insertion
    if ref_len < alt_len
        return "INS"

    # Deletion
    elseif ref_len > alt_len
        return "DEL"

    # Same length - substitution
    elseif ref_len == alt_len
        if ref_len == 1
            return "SNP"
        elseif ref_len == 2
            return "DNP"
        elseif ref_len == 3
            return "TNP"
        else
            return "ONP"
        end

    # Shouldn't reach here
    else
        return "."
    end
end

"""
    extract_hgvs_notation(transcript::TranscriptAnnotation) -> Tuple{String, String, String}

提取 HGVS 註解 (coding, protein, protein short)

# 返回
- (HGVSc, HGVSp, HGVSp_Short)
"""
function extract_hgvs_notation(transcript::TranscriptAnnotation)::Tuple{String, String, String}
    hgvsc = get(transcript, :hgvs_coding, ".")
    hgvsp = get(transcript, :hgvs_protein, ".")

    # 生成 short form (移除 p. 前綴，簡化氨基酸名稱)
    hgvsp_short = if hgvsp != "." && !isempty(hgvsp)
        simplify_hgvsp(hgvsp)
    else
        "."
    end

    return (hgvsc, hgvsp, hgvsp_short)
end

"""
    simplify_hgvsp(hgvsp::String) -> String

簡化 HGVS 蛋白質註解為短形式
例如: p.Gly12Asp → p.G12D
"""
function simplify_hgvsp(hgvsp::String)::String
    # 簡單實作：直接返回原值
    # 完整實作需要氨基酸三字母到單字母的映射
    return hgvsp
end

"""
    select_canonical_transcript(transcripts::Vector{TranscriptAnnotation}) -> Union{TranscriptAnnotation, Nothing}

根據優先級選擇 canonical transcript

# 優先級
1. RefSeq transcript (NM_ 開頭) - 通常是 canonical
2. 最嚴重的 consequence - 選擇影響最大的變異
3. 第一個 transcript (fallback)

# 參數
- `transcripts`: TranscriptAnnotation 陣列

# 返回
- 選中的 transcript，如果陣列為空則返回 nothing

# Consequence 嚴重程度 (依據 Sequence Ontology)
按從嚴重到輕微排序：
- High impact: transcript_ablation, splice_donor/acceptor, stop_gained, frameshift
- Moderate impact: inframe_indel, missense_variant
- Low impact: splice_region, synonymous_variant
- Modifier: intron_variant, upstream/downstream
"""
function select_canonical_transcript(transcripts::Vector{TranscriptAnnotation})::Union{TranscriptAnnotation, Nothing}
    if isempty(transcripts)
        return nothing
    end

    # 優先級 1: RefSeq transcript (NM_ 開頭)
    # RefSeq NM_ 是人工審核的 canonical transcript
    for trans in transcripts
        if trans.id !== nothing && startswith(trans.id, "NM_")
            return trans
        end
    end

    # 優先級 2: 根據 consequence 嚴重程度選擇
    best_trans = transcripts[1]
    best_severity = 1000  # 初始化為很大的數字

    for trans in transcripts
        if !isempty(trans.consequence)
            # 找出此 transcript 中最嚴重的 consequence
            min_severity = 1000
            for cons in trans.consequence
                severity = get(CONSEQUENCE_SEVERITY, cons, 999)  # 未知的 consequence 設為 999
                if severity < min_severity
                    min_severity = severity
                end
            end

            # 如果此 transcript 的最嚴重 consequence 比目前最好的更嚴重，則選擇它
            if min_severity < best_severity
                best_severity = min_severity
                best_trans = trans
            end
        end
    end

    return best_trans
end

"""
    extract_clinvar_info(clinvar_entries::Vector{ClinVarEntry}) -> Tuple{String, String, String, String}

提取 ClinVar 資訊

# 返回
- (ClinVar_ID, Review_Status, Significance, Disease)
"""
function extract_clinvar_info(clinvar_entries::Vector{ClinVarEntry})::Tuple{String, String, String, String}
    if isempty(clinvar_entries)
        return (".", ".", ".", ".")
    end

    # 取第一個條目
    entry = clinvar_entries[1]

    clinvar_id = entry.id !== nothing ? entry.id : "."
    review_status = entry.review_status !== nothing ? entry.review_status : "."
    significance = entry.clinical_significance !== nothing ? entry.clinical_significance : "."
    disease = isempty(entry.diseases) ? "." : join(entry.diseases, ";")

    return (clinvar_id, review_status, significance, disease)
end

"""
    extract_cosmic_id(cosmic_entries::Vector{CosmicEntry}) -> String

提取 COSMIC ID
"""
function extract_cosmic_id(cosmic_entries::Vector{CosmicEntry})::String
    if isempty(cosmic_entries)
        return "."
    end

    # 收集所有 COSMIC IDs
    ids = [e.id for e in cosmic_entries if e.id !== nothing]

    if isempty(ids)
        return "."
    end

    return join(ids, ";")
end

"""
    extract_gnomad_af(pop_freqs::Vector{PopulationFrequency}) -> Union{Float64, Nothing}

提取 gnomAD 等位基因頻率
"""
function extract_gnomad_af(pop_freqs::Vector{PopulationFrequency})::Union{Float64, Nothing}
    for pf in pop_freqs
        if pf.source == "gnomad-exome" && pf.all_af !== nothing
            return pf.all_af
        end
    end

    return nothing
end

"""
    extract_gnomad_eas_af(pop_freqs::Vector{PopulationFrequency}) -> Union{Float64, Nothing}

提取 gnomAD 東亞族群等位基因頻率
"""
function extract_gnomad_eas_af(pop_freqs::Vector{PopulationFrequency})::Union{Float64, Nothing}
    for pf in pop_freqs
        if pf.source == "gnomad-exome" && pf.eas_af !== nothing
            return pf.eas_af
        end
    end

    return nothing
end

"""
    get(transcript::TranscriptAnnotation, key::Symbol, default) -> Any

安全地從 TranscriptAnnotation 提取欄位
"""
function get(transcript::TranscriptAnnotation, key::Symbol, default)
    if key == :gene_symbol
        return transcript.gene_symbol !== nothing ? transcript.gene_symbol : default
    elseif key == :transcript_id
        return transcript.id !== nothing ? transcript.id : default
    elseif key == :consequences
        return !isempty(transcript.consequence) ? transcript.consequence : default
    elseif key == :hgvs_coding
        return transcript.hgvsc !== nothing ? transcript.hgvsc : default
    elseif key == :hgvs_protein
        return transcript.hgvsp !== nothing ? transcript.hgvsp : default
    elseif key == :is_canonical
        # TranscriptAnnotation 沒有 is_canonical 欄位，返回 default
        return default
    else
        return default
    end
end

end # module MAFConverter
