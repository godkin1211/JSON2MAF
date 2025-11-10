"""
DataStructures.jl

Define core data structures for JSON2MAF project
"""

module DataStructures

export FilterConfig, NirvanaHeader, VariantPosition, ClinVarEntry, TranscriptAnnotation
export ClinVarAssessment, QualityFilterResult, PredictiveAssessment, FilterDecision
export MAFRecord, PopulationFrequency, CosmicEntry, GeneAnnotation, NirvanaData, DataSource

# ============================================================================
# Filter configuration
# ============================================================================

"""
FilterConfig - All adjustable filtering parameters

Users can adjust these thresholds via command-line arguments
"""
Base.@kwdef struct FilterConfig
    # Quality filtering parameters
    min_total_depth::Int = 30                    # Minimum sequencing depth
    min_variant_frequency::Float64 = 0.03        # Minimum VAF (filter below this value)

    # Population frequency filtering parameters
    max_eas_af::Float64 = 0.01                   # Maximum East Asian allele frequency (1%)

    # Predictive score thresholds
    min_revel_score::Float64 = 0.75              # REVEL threshold (0.5 or 0.75)
    min_primate_ai_score::Float64 = 0.8          # PrimateAI-3D threshold
    min_dann_score::Float64 = 0.96               # DANN threshold
end

# ============================================================================
# Nirvana JSON data structures
# ============================================================================

"""
DataSource - Nirvana annotation data source
"""
struct DataSource
    name::String
    version::String
    description::Union{String, Nothing}
    release_date::Union{String, Nothing}
end

"""
NirvanaHeader - Nirvana JSON file header section
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
ClinVarEntry - ClinVar annotation entry
"""
struct ClinVarEntry
    id::Union{String, Nothing}                    # RCV ID
    allele_id::Union{String, Nothing}
    clinical_significance::Vector{String}         # Array of significance terms
    review_status::Union{String, Nothing}
    phenotypes::Vector{String}                    # Phenotypes/diseases (field name is "phenotypes" in JSON)
    last_evaluated::Union{String, Nothing}
end

"""
TranscriptAnnotation - Transcript annotation
"""
struct TranscriptAnnotation
    id::Union{String, Nothing}                    # NM_ ID
    gene_symbol::Union{String, Nothing}
    hgnc::Union{Int, Nothing}
    consequence::Vector{String}                   # SO terms
    amino_acids::Union{String, Nothing}           # e.g.: "V/M"
    cdna_pos::Union{Int, Nothing}
    cds_pos::Union{Int, Nothing}
    protein_pos::Union{Int, Nothing}
    hgvsc::Union{String, Nothing}                 # HGVS coding
    hgvsp::Union{String, Nothing}                 # HGVS protein
    is_mane_select::Union{Bool, Nothing}          # MANE Select transcript
end

"""
PopulationFrequency - Population frequency information
"""
struct PopulationFrequency
    source::String                                # "gnomad-exome", "oneKg", etc.
    all_af::Union{Float64, Nothing}               # All population frequency
    eas_af::Union{Float64, Nothing}               # East Asian population frequency
    afr_af::Union{Float64, Nothing}               # African population frequency
    amr_af::Union{Float64, Nothing}               # American population frequency
    eur_af::Union{Float64, Nothing}               # European population frequency
end

"""
CosmicEntry - COSMIC database entry
"""
struct CosmicEntry
    id::Union{String, Nothing}
    gene::Union{String, Nothing}
    mutation_type::Union{String, Nothing}
    count::Union{Int, Nothing}                    # Number of times this mutation appears in COSMIC
end

"""
VariantPosition - Complete annotation for single variant position
"""
struct VariantPosition
    chromosome::String
    start::Int
    end_pos::Int                                   # Use end_pos to avoid Julia keyword conflict
    reference_allele::String
    alternate_allele::String
    variant_type::String                           # SNV, INDEL, etc.

    # Filter status
    filters::Vector{String}                        # VCF filters field, e.g. ["PASS"] or ["filtered_reads"]

    # Sample information
    total_depth::Union{Int, Nothing}
    variant_frequencies::Union{Vector{Float64}, Nothing}

    # Annotation information
    transcripts::Vector{TranscriptAnnotation}
    clinvar::Vector{ClinVarEntry}
    cosmic::Vector{CosmicEntry}
    population_frequencies::Vector{PopulationFrequency}

    # Predictive scores (may be missing)
    primate_ai_3d::Union{Float64, Nothing}
    primate_ai::Union{Float64, Nothing}
    dann_score::Union{Float64, Nothing}
    revel_score::Union{Float64, Nothing}

    # dbSNP
    dbsnp_ids::Vector{String}
end

"""
GeneAnnotation - Gene-level annotation (auxiliary information)
"""
struct GeneAnnotation
    gene_symbol::String
    hgnc::Union{Int, Nothing}
    omim::Vector{String}
end

"""
NirvanaData - Complete Nirvana JSON parsing result
"""
struct NirvanaData
    header::NirvanaHeader
    positions::Vector{VariantPosition}
    genes::Dict{String, GeneAnnotation}
end

# ============================================================================
# Filter assessment results
# ============================================================================

"""
QualityFilterResult - Quality and population frequency filter result
"""
struct QualityFilterResult
    passes_quality::Bool
    failure_reason::Union{String, Nothing}
    depth::Union{Int, Nothing}
    variant_frequency::Union{Float64, Nothing}
    eas_allele_frequency::Union{Float64, Nothing}
end

"""
ClinVarAssessment - ClinVar pathogenicity assessment result
"""
struct ClinVarAssessment
    is_pathogenic::Bool
    is_likely_pathogenic::Bool
    selected_entry::Union{ClinVarEntry, Nothing}
    confidence_level::String                       # "high", "medium", "low"
    reason::String
end

"""
PredictiveAssessment - Predictive score assessment result
"""
struct PredictiveAssessment
    suggests_pathogenic::Bool
    contributing_scores::Dict{String, Float64}     # Scores exceeding threshold
    confidence::Float64
    support_count::Int                             # How many scores support
    has_primate_ai_support::Bool                   # Whether PrimateAI-3D supports
end

"""
FilterDecision - Final filtering decision
"""
struct FilterDecision
    should_include::Bool
    pathogenicity_class::String                    # "Pathogenic", "Likely pathogenic", "Excluded"
    primary_evidence::String                       # "ClinVar", "Predictive", "None"
    justification::String                          # Detailed explanation of decision rationale
end

# ============================================================================
# MAF format
# ============================================================================

"""
MAFRecord - Single variant record in MAF format
"""
Base.@kwdef struct MAFRecord
    # Basic information
    hugo_symbol::String = ""
    chromosome::String = ""
    start_position::Int = 0
    end_position::Int = 0
    strand::String = "+"

    # Variant classification
    variant_classification::String = ""
    variant_type::String = ""

    # Alleles
    reference_allele::String = ""
    tumor_seq_allele1::String = ""                 # Usually same as reference
    tumor_seq_allele2::String = ""                 # Variant allele

    # Sample information
    tumor_sample_barcode::String = ""

    # HGVS naming
    hgvsc::String = ""
    hgvsp::String = ""
    hgvsp_short::String = ""

    # Transcript
    transcript_id::String = ""

    # Database IDs
    dbsnp_rs::String = ""
    dbsnp_val_status::String = ""
    cosmic_id::String = ""

    # ClinVar information
    clinvar_id::String = ""
    clinvar_review_status::String = ""
    clinvar_significance::String = ""
    clinvar_disease::String = ""

    # Predictive scores (use nothing to indicate missing)
    primate_ai_score::Union{String, Nothing} = nothing
    dann_score::Union{String, Nothing} = nothing
    revel_score::Union{String, Nothing} = nothing

    # Population frequency (use nothing to indicate missing)
    gnomad_af::Union{String, Nothing} = nothing
    gnomad_eas_af::Union{String, Nothing} = nothing

    # Sequencing quality
    depth::String = ""
    vaf::String = ""
end

end # module DataStructures
