"""
MAFConverter.jl

MAF format conversion module
Converts Nirvana annotated variants to MAF (Mutation Annotation Format) format
"""

module MAFConverter

using ..DataStructures
import Base: get

export variant_to_maf, map_variant_classification, map_variant_type,
       extract_hgvs_notation, select_canonical_transcript

# Global constants: avoid repeated allocation of empty strings
const EMPTY_STRING = ""
const DOT_STRING = "."

# Consequence severity mapping (global constant to avoid recreating Dict on each function call)
# Lower number means more severe, based on Sequence Ontology and VEP classification
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

Convert variant position to MAF record

# Parameters
- `variant`: Variant position data
- `decision`: Filter decision result

# Returns
- `MAFRecord`: MAF format record

# Example
```julia
maf_record = variant_to_maf(variant, decision)
```
"""
function variant_to_maf(variant::VariantPosition, decision::FilterDecision)::MAFRecord
    # Select canonical transcript
    canonical = select_canonical_transcript(variant.transcripts)

    # Extract gene symbol
    hugo_symbol = canonical !== nothing ? get(canonical, :gene_symbol, ".") : "."

    # Extract HGVS
    hgvsc, hgvsp, hgvsp_short = if canonical !== nothing
        extract_hgvs_notation(canonical)
    else
        (".", ".", ".")
    end

    # Extract Transcript ID
    transcript_id = canonical !== nothing ? get(canonical, :transcript_id, ".") : "."

    # Map variant classification
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

    # Map variant type
    variant_type_maf = map_variant_type(variant.reference_allele, variant.alternate_allele)

    # Extract ClinVar information
    clinvar_id, clinvar_review, clinvar_sig, clinvar_disease = extract_clinvar_info(variant.clinvar)

    # Extract COSMIC ID
    cosmic_id = extract_cosmic_id(variant.cosmic)

    # Extract dbSNP
    dbsnp_rs = isempty(variant.dbsnp_ids) ? "." : variant.dbsnp_ids[1]

    # Extract gnomAD AF
    gnomad_af = extract_gnomad_af(variant.population_frequencies)

    # Extract predictive scores and convert to string (use nothing for missing values)
    # PrimateAI: prioritize 3D version
    primate_ai = variant.primate_ai_3d !== nothing ? variant.primate_ai_3d :
                 (variant.primate_ai !== nothing ? variant.primate_ai : nothing)
    primate_ai_str = isnothing(primate_ai) ? nothing : string(primate_ai)
    dann_str = isnothing(variant.dann_score) ? nothing : string(variant.dann_score)
    revel_str = isnothing(variant.revel_score) ? nothing : string(variant.revel_score)
    gnomad_af_str = isnothing(gnomad_af) ? nothing : string(gnomad_af)

    # Extract East Asian AF
    gnomad_eas_af = extract_gnomad_eas_af(variant.population_frequencies)
    gnomad_eas_af_str = isnothing(gnomad_eas_af) ? nothing : string(gnomad_eas_af)

    # Extract depth and VAF - use global constants to reduce allocations
    depth_str = isnothing(variant.total_depth) ? EMPTY_STRING : string(variant.total_depth)
    vaf_str = (isnothing(variant.variant_frequencies) || isempty(variant.variant_frequencies)) ?
              EMPTY_STRING : string(variant.variant_frequencies[1])

    return MAFRecord(
        hugo_symbol = hugo_symbol,
        chromosome = variant.chromosome,
        start_position = variant.start,
        end_position = variant.end_pos,
        strand = "+",  # Nirvana usually doesn't provide, default to +
        variant_classification = variant_classification,
        variant_type = variant_type_maf,
        reference_allele = variant.reference_allele,
        tumor_seq_allele1 = variant.reference_allele,  # Usually same as reference
        tumor_seq_allele2 = variant.alternate_allele,  # Variant allele
        tumor_sample_barcode = "",  # Need to extract from samples, currently empty
        hgvsc = hgvsc,
        hgvsp = hgvsp,
        hgvsp_short = hgvsp_short,
        transcript_id = transcript_id,
        dbsnp_rs = dbsnp_rs,
        dbsnp_val_status = "",  # Not available
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

Map Nirvana consequence to MAF Variant_Classification

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
        return consequence  # Return original value
    end
end

"""
    map_variant_type(ref::String, alt::String) -> String

Map variant type to MAF Variant_Type

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

Extract HGVS notation (coding, protein, protein short)

# Returns
- (HGVSc, HGVSp, HGVSp_Short)
"""
function extract_hgvs_notation(transcript::TranscriptAnnotation)::Tuple{String, String, String}
    hgvsc = get(transcript, :hgvs_coding, ".")
    hgvsp = get(transcript, :hgvs_protein, ".")

    # Generate short form (remove p. prefix, simplify amino acid names)
    hgvsp_short = if hgvsp != "." && !isempty(hgvsp)
        simplify_hgvsp(hgvsp)
    else
        "."
    end

    return (hgvsc, hgvsp, hgvsp_short)
end

"""
    simplify_hgvsp(hgvsp::String) -> String

Simplify HGVS protein notation to short form
Example: p.Gly12Asp → p.G12D
"""
const AMINO_ACID_MAP = Dict(
    "Ala" => "A", "Arg" => "R", "Asn" => "N", "Asp" => "D", "Cys" => "C",
    "Glu" => "E", "Gln" => "Q", "Gly" => "G", "His" => "H", "Ile" => "I",
    "Leu" => "L", "Lys" => "K", "Met" => "M", "Phe" => "F", "Pro" => "P",
    "Ser" => "S", "Thr" => "T", "Trp" => "W", "Tyr" => "Y", "Val" => "V",
    "Ter" => "*"
)

"""
    simplify_hgvsp(hgvsp::String) -> String

Simplify HGVS protein notation to short form
Example: p.Gly12Asp → p.G12D
"""
function simplify_hgvsp(hgvsp::String)::String
    if !startswith(hgvsp, "p.")
        return hgvsp
    end

    # Substitution: p.Gly12Asp
    sub_match = match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", hgvsp)
    if sub_match !== nothing
        ref = get(AMINO_ACID_MAP, sub_match.captures[1], sub_match.captures[1])
        pos = sub_match.captures[2]
        alt = get(AMINO_ACID_MAP, sub_match.captures[3], sub_match.captures[3])
        return "p.$(ref)$(pos)$(alt)"
    end

    # Deletion: p.Gly12del
    del_match = match(r"p\.([A-Z][a-z]{2})(\d+)del", hgvsp)
    if del_match !== nothing
        ref = get(AMINO_ACID_MAP, del_match.captures[1], del_match.captures[1])
        pos = del_match.captures[2]
        return "p.$(ref)$(pos)del"
    end

    # Insertion: p.Gly12_Ala13insAsp
    ins_match = match(r"p\.([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)ins([A-Z][a-z]{2})", hgvsp)
    if ins_match !== nothing
        ref1 = get(AMINO_ACID_MAP, ins_match.captures[1], ins_match.captures[1])
        pos1 = ins_match.captures[2]
        ref2 = get(AMINO_ACID_MAP, ins_match.captures[3], ins_match.captures[3])
        pos2 = ins_match.captures[4]
        ins = get(AMINO_ACID_MAP, ins_match.captures[5], ins_match.captures[5])
        return "p.$(ref1)$(pos1)_$(ref2)$(pos2)ins$(ins)"
    end

    return hgvsp  # Return original if no match
end

"""
    select_canonical_transcript(transcripts::Vector{TranscriptAnnotation}) -> Union{TranscriptAnnotation, Nothing}

Select canonical transcript based on priority

# Priority
1. RefSeq transcript (starting with NM_) - usually canonical
2. Most severe consequence - choose the most impactful variant
3. First transcript (fallback)

# Parameters
- `transcripts`: Array of TranscriptAnnotation

# Returns
- Selected transcript, or nothing if array is empty

# Consequence severity (based on Sequence Ontology)
Sorted from severe to mild:
- High impact: transcript_ablation, splice_donor/acceptor, stop_gained, frameshift
- Moderate impact: inframe_indel, missense_variant
- Low impact: splice_region, synonymous_variant
- Modifier: intron_variant, upstream/downstream
"""
function select_canonical_transcript(transcripts::Vector{TranscriptAnnotation})::Union{TranscriptAnnotation, Nothing}
    if isempty(transcripts)
        return nothing
    end

    # Priority 1: MANE Select transcript
    mane_transcripts = filter(t -> get(t, :is_mane_select, false), transcripts)
    if !isempty(mane_transcripts)
        return mane_transcripts[1]
    end

    # Priority 2: RefSeq transcript (starting with NM_)
    # RefSeq NM_ is manually reviewed canonical transcript
    for trans in transcripts
        if trans.id !== nothing && startswith(trans.id, "NM_")
            return trans
        end
    end

    # Priority 3: Select based on consequence severity
    best_trans = transcripts[1]
    best_severity = 1000  # Initialize with a large number

    for trans in transcripts
        if !isempty(trans.consequence)
            # Find the most severe consequence in this transcript
            min_severity = 1000
            for cons in trans.consequence
                severity = get(CONSEQUENCE_SEVERITY, cons, 999)  # Unknown consequences set to 999
                if severity < min_severity
                    min_severity = severity
                end
            end

            # If this transcript's most severe consequence is more severe than the current best, select it
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

Extract ClinVar information

# Returns
- (ClinVar_ID, Review_Status, Significance, Disease)
"""
function extract_clinvar_info(clinvar_entries::Vector{ClinVarEntry})::Tuple{String, String, String, String}
    if isempty(clinvar_entries)
        return (".", ".", ".", ".")
    end

    # Take the first entry
    entry = clinvar_entries[1]

    clinvar_id = entry.id !== nothing ? entry.id : "."
    review_status = entry.review_status !== nothing ? entry.review_status : "."
    # Join array elements with ", "
    significance = !isempty(entry.clinical_significance) ? join(entry.clinical_significance, ", ") : "."
    disease = isempty(entry.phenotypes) ? "." : join(entry.phenotypes, ";")

    return (clinvar_id, review_status, significance, disease)
end

"""
    extract_cosmic_id(cosmic_entries::Vector{CosmicEntry}) -> String

Extract COSMIC ID
"""
function extract_cosmic_id(cosmic_entries::Vector{CosmicEntry})::String
    if isempty(cosmic_entries)
        return "."
    end

    # Collect all COSMIC IDs
    ids = [e.id for e in cosmic_entries if e.id !== nothing]

    if isempty(ids)
        return "."
    end

    return join(ids, ";")
end

"""
    extract_gnomad_af(pop_freqs::Vector{PopulationFrequency}) -> Union{Float64, Nothing}

Extract gnomAD allele frequency
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

Extract gnomAD East Asian population allele frequency
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

Safely extract fields from TranscriptAnnotation
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
    elseif key == :is_mane_select
        return transcript.is_mane_select !== nothing ? transcript.is_mane_select : default
    elseif key == :is_canonical
        # TranscriptAnnotation doesn't have is_canonical field, return default
        return default
    else
        return default
    end
end

end # module MAFConverter
