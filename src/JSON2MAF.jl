"""
JSON2MAF.jl

Main module for Nirvana JSON pathogenic variant filtering tool

Filters Illumina Nirvana annotation results to retain Pathogenic/Likely Pathogenic variants,
and outputs in MAF (Mutation Annotation Format) format.
"""

module JSON2MAF

# Version information
const VERSION = v"0.4.0"

# Export core modules
include("utils/DataStructures.jl")
using .DataStructures

include("utils/FilterConfig.jl")
using .FilterConfigModule

# Parser module
include("parser/NirvanaParser.jl")
using .NirvanaParser

# Filter modules
include("filters/QualityFilter.jl")
using .QualityFilter

include("filters/ClinVarFilter.jl")
using .ClinVarFilter

include("filters/PredictiveScores.jl")
using .PredictiveScores

include("filters/DecisionEngine.jl")
using .DecisionEngine

# Converter module
include("converters/MAFConverter.jl")
using .MAFConverter

# Writer module
include("writers/MAFWriter.jl")
using .MAFWriter

# Export all public APIs
export FilterConfig, create_filter_config, validate_config, display_config
export NirvanaHeader, VariantPosition, ClinVarEntry, TranscriptAnnotation
export ClinVarAssessment, QualityFilterResult, PredictiveAssessment, FilterDecision
export MAFRecord, PopulationFrequency, CosmicEntry, GeneAnnotation, NirvanaData, DataSource
export process_nirvana_parallel_with_stats
export apply_quality_filters, check_sequencing_quality, check_population_frequency
export assess_clinvar_pathogenicity, get_review_status_priority, is_cancer_related, resolve_conflicting_entries
export assess_predictive_scores, get_primate_ai_score, get_dann_score, get_revel_score, is_in_cosmic
export make_filter_decision, count_supporting_predictive_scores, has_primate_ai_support
export variant_to_maf, map_variant_classification, map_variant_type, extract_hgvs_notation, select_canonical_transcript
export maf_header, mafrecord_to_row
export create_maf_writer, write_maf_batch, close_maf_writer, get_total_written
export merge_maf_files, count_maf_records

"""
    version()

Returns the JSON2MAF version number
"""
version() = VERSION

end # module JSON2MAF
