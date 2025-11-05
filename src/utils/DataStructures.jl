"""
DataStructures.jl

定義 JSON2MAF 專案的核心資料結構
"""

module DataStructures

export FilterConfig, NirvanaHeader, VariantPosition, ClinVarEntry, TranscriptAnnotation
export ClinVarAssessment, QualityFilterResult, PredictiveAssessment, FilterDecision
export MAFRecord, PopulationFrequency, CosmicEntry, GeneAnnotation, NirvanaData, DataSource

# ============================================================================
# 過濾配置
# ============================================================================

"""
FilterConfig - 所有可調整的過濾參數

使用者可以通過命令行參數調整這些閾值
"""
Base.@kwdef struct FilterConfig
    # 品質過濾參數
    min_total_depth::Int = 30                    # 最小測序深度
    min_variant_frequency::Float64 = 0.03        # 最小 VAF (過濾低於此值)

    # 族群頻率過濾參數
    max_eas_af::Float64 = 0.01                   # 東亞族群最大等位基因頻率 (1%)

    # 預測分數閾值
    min_revel_score::Float64 = 0.75              # REVEL閾值 (0.5或0.75)
    min_primate_ai_score::Float64 = 0.8          # PrimateAI-3D閾值
    min_dann_score::Float64 = 0.96               # DANN閾值
end

# ============================================================================
# Nirvana JSON 資料結構
# ============================================================================

"""
DataSource - Nirvana 註解資料來源
"""
struct DataSource
    name::String
    version::String
    description::Union{String, Nothing}
    release_date::Union{String, Nothing}
end

"""
NirvanaHeader - Nirvana JSON 檔案的 header 區段
"""
struct NirvanaHeader
    annotator::String
    creation_time::String
    genome_assembly::String
    schema_version::Int
    data_sources::Vector{DataSource}
    samples::Vector{String}
end

"""
ClinVarEntry - ClinVar 註解條目
"""
struct ClinVarEntry
    id::Union{String, Nothing}                    # RCV ID
    allele_id::Union{String, Nothing}
    clinical_significance::Union{String, Nothing}
    review_status::Union{String, Nothing}
    diseases::Vector{String}
    last_evaluated::Union{String, Nothing}
end

"""
TranscriptAnnotation - 轉錄本註解
"""
struct TranscriptAnnotation
    id::Union{String, Nothing}                    # NM_ ID
    gene_symbol::Union{String, Nothing}
    hgnc::Union{Int, Nothing}
    consequence::Vector{String}                   # SO terms
    amino_acids::Union{String, Nothing}           # 例如: "V/M"
    cdna_pos::Union{Int, Nothing}
    cds_pos::Union{Int, Nothing}
    protein_pos::Union{Int, Nothing}
    hgvsc::Union{String, Nothing}                 # HGVS coding
    hgvsp::Union{String, Nothing}                 # HGVS protein
end

"""
PopulationFrequency - 族群頻率資訊
"""
struct PopulationFrequency
    source::String                                # "gnomad-exome", "oneKg", etc.
    all_af::Union{Float64, Nothing}               # 全體族群頻率
    eas_af::Union{Float64, Nothing}               # 東亞族群頻率
    afr_af::Union{Float64, Nothing}               # 非洲族群頻率
    amr_af::Union{Float64, Nothing}               # 美洲族群頻率
    eur_af::Union{Float64, Nothing}               # 歐洲族群頻率
end

"""
CosmicEntry - COSMIC 資料庫條目
"""
struct CosmicEntry
    id::Union{String, Nothing}
    gene::Union{String, Nothing}
    mutation_type::Union{String, Nothing}
    count::Union{Int, Nothing}                    # 該突變在COSMIC中的出現次數
end

"""
VariantPosition - 單一變異位點的完整註解
"""
struct VariantPosition
    chromosome::String
    start::Int
    end_pos::Int                                   # 使用 end_pos 避免與 Julia 關鍵字衝突
    reference_allele::String
    alternate_allele::String
    variant_type::String                           # SNV, INDEL, etc.

    # 過濾狀態
    filters::Vector{String}                        # VCF filters 欄位，例如 ["PASS"] 或 ["filtered_reads"]

    # 樣本資訊
    total_depth::Union{Int, Nothing}
    variant_frequencies::Union{Vector{Float64}, Nothing}

    # 註解資訊
    transcripts::Vector{TranscriptAnnotation}
    clinvar::Vector{ClinVarEntry}
    cosmic::Vector{CosmicEntry}
    population_frequencies::Vector{PopulationFrequency}

    # 預測分數 (可能缺失)
    primate_ai_3d::Union{Float64, Nothing}
    primate_ai::Union{Float64, Nothing}
    dann_score::Union{Float64, Nothing}
    revel_score::Union{Float64, Nothing}

    # dbSNP
    dbsnp_ids::Vector{String}
end

"""
GeneAnnotation - 基因層級的註解 (輔助資訊)
"""
struct GeneAnnotation
    gene_symbol::String
    hgnc::Union{Int, Nothing}
    omim::Vector{String}
end

"""
NirvanaData - 完整的 Nirvana JSON 解析結果
"""
struct NirvanaData
    header::NirvanaHeader
    positions::Vector{VariantPosition}
    genes::Dict{String, GeneAnnotation}
end

# ============================================================================
# 過濾評估結果
# ============================================================================

"""
QualityFilterResult - 品質與族群頻率過濾結果
"""
struct QualityFilterResult
    passes_quality::Bool
    failure_reason::Union{String, Nothing}
    depth::Union{Int, Nothing}
    variant_frequency::Union{Float64, Nothing}
    eas_allele_frequency::Union{Float64, Nothing}
end

"""
ClinVarAssessment - ClinVar 致病性評估結果
"""
struct ClinVarAssessment
    is_pathogenic::Bool
    is_likely_pathogenic::Bool
    selected_entry::Union{ClinVarEntry, Nothing}
    confidence_level::String                       # "high", "medium", "low"
    reason::String
end

"""
PredictiveAssessment - 預測分數評估結果
"""
struct PredictiveAssessment
    suggests_pathogenic::Bool
    contributing_scores::Dict{String, Float64}     # 超過閾值的分數
    confidence::Float64
    support_count::Int                             # 有多少個分數支持
    has_primate_ai_support::Bool                   # PrimateAI-3D 是否支持
end

"""
FilterDecision - 最終過濾決策
"""
struct FilterDecision
    should_include::Bool
    pathogenicity_class::String                    # "Pathogenic", "Likely pathogenic", "Excluded"
    primary_evidence::String                       # "ClinVar", "Predictive", "None"
    justification::String                          # 詳細說明決策理由
end

# ============================================================================
# MAF 格式
# ============================================================================

"""
MAFRecord - MAF 格式的單一變異記錄
"""
Base.@kwdef struct MAFRecord
    # 基本資訊
    hugo_symbol::String = ""
    chromosome::String = ""
    start_position::Int = 0
    end_position::Int = 0
    strand::String = "+"

    # 變異分類
    variant_classification::String = ""
    variant_type::String = ""

    # 等位基因
    reference_allele::String = ""
    tumor_seq_allele1::String = ""                 # 通常與 reference 相同
    tumor_seq_allele2::String = ""                 # 變異等位基因

    # 樣本資訊
    tumor_sample_barcode::String = ""

    # HGVS 命名
    hgvsc::String = ""
    hgvsp::String = ""
    hgvsp_short::String = ""

    # 轉錄本
    transcript_id::String = ""

    # 資料庫 ID
    dbsnp_rs::String = ""
    dbsnp_val_status::String = ""
    cosmic_id::String = ""

    # ClinVar 資訊
    clinvar_id::String = ""
    clinvar_review_status::String = ""
    clinvar_significance::String = ""
    clinvar_disease::String = ""

    # 預測分數 (使用 nothing 表示缺失)
    primate_ai_score::Union{String, Nothing} = nothing
    dann_score::Union{String, Nothing} = nothing
    revel_score::Union{String, Nothing} = nothing

    # 族群頻率 (使用 nothing 表示缺失)
    gnomad_af::Union{String, Nothing} = nothing
    gnomad_eas_af::Union{String, Nothing} = nothing

    # 測序品質
    depth::String = ""
    vaf::String = ""
end

end # module DataStructures
