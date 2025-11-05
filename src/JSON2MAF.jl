"""
JSON2MAF.jl

Nirvana JSON 致癌性重要位點過濾工具主模組

將 Illumina Nirvana 註解結果過濾，保留 Pathogenic/Likely Pathogenic 變異，
並輸出為 MAF (Mutation Annotation Format) 格式。
"""

module JSON2MAF

# 版本資訊
const VERSION = v"0.4.0"

# 匯出核心模組
include("utils/DataStructures.jl")
using .DataStructures

include("utils/FilterConfig.jl")
using .FilterConfigModule

# 解析模組
include("parser/NirvanaParser.jl")
using .NirvanaParser

# 過濾模組
include("filters/QualityFilter.jl")
using .QualityFilter

include("filters/ClinVarFilter.jl")
using .ClinVarFilter

include("filters/PredictiveScores.jl")
using .PredictiveScores

include("filters/DecisionEngine.jl")
using .DecisionEngine

# 轉換模組
include("converters/MAFConverter.jl")
using .MAFConverter

# 寫入模組
include("writers/MAFWriter.jl")
using .MAFWriter

# 匯出所有公開 API
export FilterConfig, create_filter_config, validate_config, display_config
export NirvanaHeader, VariantPosition, ClinVarEntry, TranscriptAnnotation
export ClinVarAssessment, QualityFilterResult, PredictiveAssessment, FilterDecision
export MAFRecord, PopulationFrequency, CosmicEntry, GeneAnnotation, NirvanaData, DataSource
export parse_nirvana_json, parse_nirvana_with_prefilter, process_nirvana_streaming, process_nirvana_parallel, process_nirvana_parallel_no_prefilter, process_nirvana_parallel_with_stats, process_nirvana_streaming_parallel
export apply_quality_filters, check_sequencing_quality, check_population_frequency
export assess_clinvar_pathogenicity, get_review_status_priority, is_cancer_related, resolve_conflicting_entries
export assess_predictive_scores, get_primate_ai_score, get_dann_score, get_revel_score, is_in_cosmic
export make_filter_decision, count_supporting_predictive_scores, has_primate_ai_support
export variant_to_maf, map_variant_classification, map_variant_type, extract_hgvs_notation, select_canonical_transcript
export write_maf_file, maf_header, mafrecord_to_row
export create_maf_writer, write_maf_batch, close_maf_writer, get_total_written
export merge_maf_files, count_maf_records

"""
    version()

返回 JSON2MAF 版本號
"""
version() = VERSION

end # module JSON2MAF
