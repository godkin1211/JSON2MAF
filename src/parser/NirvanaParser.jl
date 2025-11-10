"""
NirvanaParser.jl

Nirvana JSON file parsing module
Supports efficient parsing of gzipped JSON files
"""

module NirvanaParser

using JSON3
using CodecZlib
using ProgressMeter
using ..DataStructures
using ..FilterConfigModule

export process_nirvana_parallel_with_stats

"""
    parse_header(header_json) -> NirvanaHeader

Parse Nirvana JSON header section
"""
function parse_header(header_json)::NirvanaHeader
    data_sources = [
        DataSource(
            ds.name,
            ds.version,
            get(ds, :description, nothing),
            get(ds, :releaseDate, nothing)
        )
        for ds in header_json.dataSources
    ]

    return NirvanaHeader(
        header_json.annotator,
        header_json.creationTime,
        header_json.genomeAssembly,
        header_json.schemaVersion,
        data_sources,
        collect(header_json.samples)
    )
end

"""
    parse_position(pos_json) -> Union{VariantPosition, Nothing}

Parse a single variant position
"""
function parse_position(pos_json)::Union{VariantPosition, Nothing}
    # Basic information
    chromosome = pos_json.chromosome
    position = pos_json.position
    ref_allele = pos_json.refAllele
    
    if !haskey(pos_json, :altAlleles) || isempty(pos_json.altAlleles)
        @warn "No altAlleles found for position $(chromosome):$(position). Skipping."
        return nothing
    end
    alt_alleles = collect(pos_json.altAlleles)

    # Parse filters field
    filters = haskey(pos_json, :filters) ? collect(String, pos_json.filters) : String["PASS"]

    # Sample information - take first sample's data
    if haskey(pos_json, :samples) && length(pos_json.samples) > 1
        @warn "Multiple samples found at position $(chromosome):$(position). Only the first sample will be processed."
    end
    sample = if haskey(pos_json, :samples) && !isempty(pos_json.samples)
        pos_json.samples[1]
    else
        nothing
    end

    total_depth = sample !== nothing && haskey(sample, :totalDepth) ?
                  sample.totalDepth : nothing
    variant_frequencies = sample !== nothing && haskey(sample, :variantFrequencies) ?
                         collect(Float64, sample.variantFrequencies) : nothing

    # Variant annotations - process first variant (usually only one altAllele)
    if !haskey(pos_json, :variants) || isempty(pos_json.variants)
        @warn "No variant annotations found for position $(chromosome):$(position). Skipping."
        return nothing
    end
    variants_data = pos_json.variants
    variant = variants_data[1]

    # Parse transcripts
    transcripts = haskey(variant, :transcripts) ?
                 parse_transcripts(variant.transcripts) :
                 TranscriptAnnotation[]

    # Parse ClinVar
    clinvar = haskey(variant, :clinvar) ?
             parse_clinvar_entries(variant.clinvar) :
             ClinVarEntry[]

    # Parse COSMIC
    cosmic = haskey(variant, :cosmic) ?
            parse_cosmic_entries(variant.cosmic) :
            CosmicEntry[]

    # Parse population frequencies
    pop_freqs = parse_population_frequencies(variant)

    # Extract predictive scores
    primate_ai_3d = extract_primate_ai_3d_score(variant)
    primate_ai = extract_primate_ai_score(variant)
    dann_score = haskey(variant, :dannScore) ? Float64(variant.dannScore) : nothing
    revel_score = extract_revel_score(variant)

    # dbSNP IDs
    dbsnp_ids = haskey(variant, :dbsnp) ? collect(String, variant.dbsnp) : String[]

    return VariantPosition(
        chromosome,
        position,
        position,  # end position (same for SNVs)
        ref_allele,
        alt_alleles[1],  # First alternate allele
        variant.variantType,
        filters,  # VCF filters field
        total_depth,
        variant_frequencies,
        transcripts,
        clinvar,
        cosmic,
        pop_freqs,
        primate_ai_3d,
        primate_ai,
        dann_score,
        revel_score,
        dbsnp_ids
    )
end

"""
    parse_position_field(val) -> Union{Int, Nothing}

Parse position field, which may be an integer, numeric string, or string containing a range (e.g., "303/2618")
Takes the first number
"""
function parse_position_field(val)::Union{Int, Nothing}
    if val === nothing
        return nothing
    elseif val isa Integer
        return val
    elseif val isa String
        # Handle complex position formats like "1535-1539/4210"
        # First extract the part before "/" if present
        working_val = val
        if contains(val, "/")
            working_val = split(val, "/")[1]
        end

        # Then extract the first number if it's a range with "-"
        if contains(working_val, "-")
            working_val = split(working_val, "-")[1]
        end

        # Finally try to parse the cleaned value
        parsed_val = tryparse(Int, working_val)

        if parsed_val === nothing
            @warn "Could not parse position field: $(val)"
        end
        return parsed_val
    else
        return nothing
    end
end

"""
    parse_transcripts(transcripts_json) -> Vector{TranscriptAnnotation}

Parse transcript annotations
"""
function parse_transcripts(transcripts_json)::Vector{TranscriptAnnotation}
    result = TranscriptAnnotation[]

    for t in transcripts_json
        consequence = haskey(t, :consequence) ? collect(String, t.consequence) : String[]

        # hgnc may be an integer ID or string (gene symbol)
        hgnc_val = get(t, :hgnc, nothing)
        hgnc_id = if hgnc_val !== nothing
            if hgnc_val isa Integer
                hgnc_val
            elseif hgnc_val isa String
                # Try to parse as integer, return nothing if failed
                tryparse(Int, hgnc_val)
            else
                nothing
            end
        else
            nothing
        end

        # Parse position fields (may contain ranges or multiple values)
        cdna_pos = parse_position_field(get(t, :cdnaPos, nothing))
        cds_pos = parse_position_field(get(t, :cdsPos, nothing))
        protein_pos = parse_position_field(get(t, :proteinPos, nothing))

        push!(result, TranscriptAnnotation(
            get(t, :transcript, nothing),
            get(t, :hgnc, nothing) isa String ? get(t, :hgnc, nothing) : nothing,
            hgnc_id,
            consequence,
            get(t, :aminoAcids, nothing),
            cdna_pos,
            cds_pos,
            protein_pos,
            get(t, :hgvsc, nothing),
            get(t, :hgvsp, nothing),
            get(t, :isManeSelect, nothing)
        ))
    end

    return result
end

"""
    parse_clinvar_entries(clinvar_json) -> Vector{ClinVarEntry}

Parse ClinVar annotations
"""
function parse_clinvar_entries(clinvar_json)::Vector{ClinVarEntry}
    result = ClinVarEntry[]

    for cv in clinvar_json
        # Parse phenotypes array (field name is "phenotypes", not "diseases")
        phenotypes = haskey(cv, :phenotypes) ? collect(String, cv.phenotypes) : String[]

        # Parse significance as array (field name is "significance", not "clinicalSignificance")
        significance = haskey(cv, :significance) ? collect(String, cv.significance) : String[]

        push!(result, ClinVarEntry(
            get(cv, :id, nothing),
            get(cv, :alleleId, nothing) !== nothing ? string(get(cv, :alleleId, "")) : nothing,
            significance,
            get(cv, :reviewStatus, nothing),
            phenotypes,
            get(cv, :lastEvaluated, nothing)
        ))
    end

    return result
end

"""
    parse_cosmic_entries(cosmic_json) -> Vector{CosmicEntry}

Parse COSMIC annotations
"""
function parse_cosmic_entries(cosmic_json)::Vector{CosmicEntry}
    result = CosmicEntry[]

    for c in cosmic_json
        push!(result, CosmicEntry(
            get(c, :id, nothing),
            get(c, :gene, nothing),
            get(c, :mutationType, nothing),
            get(c, :numSamples, nothing)
        ))
    end

    return result
end

"""
    parse_population_frequencies(variant_json) -> Vector{PopulationFrequency}

Extract population frequency information from variant annotations
"""
function parse_population_frequencies(variant_json)::Vector{PopulationFrequency}
    result = PopulationFrequency[]

    # gnomAD
    if haskey(variant_json, :gnomad)
        gn = variant_json.gnomad
        push!(result, PopulationFrequency(
            "gnomad",
            get(gn, :allAf, nothing),
            get(gn, :easAf, nothing),
            get(gn, :afrAf, nothing),
            get(gn, :amrAf, nothing),
            get(gn, :nfeAf, nothing)
        ))
    end

    # gnomAD-exome
    if haskey(variant_json, Symbol("gnomad-exome"))
        gne = variant_json[Symbol("gnomad-exome")]
        push!(result, PopulationFrequency(
            "gnomad-exome",
            get(gne, :allAf, nothing),
            get(gne, :easAf, nothing),
            get(gne, :afrAf, nothing),
            get(gne, :amrAf, nothing),
            get(gne, :nfeAf, nothing)
        ))
    end

    # 1000 Genomes (oneKg)
    if haskey(variant_json, :oneKg)
        okg = variant_json.oneKg
        push!(result, PopulationFrequency(
            "oneKg",
            get(okg, :allAf, nothing),
            get(okg, :easAf, nothing),
            get(okg, :afrAf, nothing),
            get(okg, :amrAf, nothing),
            get(okg, :eurAf, nothing)
        ))
    end

    return result
end

"""
    extract_primate_ai_3d_score(variant_json) -> Union{Float64, Nothing}

Extract PrimateAI-3D score
"""
function extract_primate_ai_3d_score(variant_json)::Union{Float64, Nothing}
    if haskey(variant_json, Symbol("primateAI-3D"))
        entries = variant_json[Symbol("primateAI-3D")]
        if length(entries) > 0 && haskey(entries[1], :score)
            return Float64(entries[1].score)
        end
    end
    return nothing
end

"""
    extract_primate_ai_score(variant_json) -> Union{Float64, Nothing}

Extract PrimateAI score
"""
function extract_primate_ai_score(variant_json)::Union{Float64, Nothing}
    if haskey(variant_json, :primateAI)
        entries = variant_json.primateAI
        if length(entries) > 0 && haskey(entries[1], :scorePercentile)
            return Float64(entries[1].scorePercentile)
        end
    end
    return nothing
end

"""
    extract_revel_score(variant_json) -> Union{Float64, Nothing}

Extract REVEL score
"""
function extract_revel_score(variant_json)::Union{Float64, Nothing}
    if haskey(variant_json, :revel) && haskey(variant_json.revel, :score)
        return Float64(variant_json.revel.score)
    end
    return nothing
end

"""
    process_nirvana_parallel_with_stats(callback, filepath, config, thread_stats; num_threads)

Multi-threaded parallel processing of Nirvana JSON file (prefilter version with statistics tracking)

This version tracks statistics during the prefiltering stage, solving the incomplete statistics problem
of the original `process_nirvana_parallel`, while maintaining prefilter performance advantages.

# Parameters
- `callback`: Processing function `(variant, idx, thread_id) -> Nothing`
- `filepath`: Path to Nirvana JSON.gz file
- `config`: Filter configuration parameters
- `thread_stats`: Array of statistics dictionaries for each thread (will be updated with prefilter statistics in this function)
- `num_threads`: Number of threads (default uses Threads.nthreads())

# thread_stats format
Each thread needs a Dict with at least the following keys:
- `:failed_depth` - Number of variants with insufficient depth
- `:failed_vaf` - Number of variants with VAF too low
- `:failed_af` - Number of variants with population frequency too high

# Returns
- `NirvanaHeader`: Only returns header

# Features
- Records filtered variant statistics during prefiltering stage
- Maintains prefilter performance advantages
- Provides complete statistical tracking
"""
function process_nirvana_parallel_with_stats(
    callback::Function,
    filepath::String,
    config::FilterConfig,
    thread_stats::Vector{Dict{Symbol, Int}};
    num_threads::Union{Int, Nothing}=nothing
)::NirvanaHeader

    if !isfile(filepath)
        error("File not found: $filepath")
    end

    # Get number of threads
    nthreads = num_threads === nothing ? Threads.nthreads() : num_threads

    if nthreads == 1
        @warn "Only 1 thread available. Consider setting JULIA_NUM_THREADS for better performance."
    end

    # Read and parse JSON
    json_data = if endswith(filepath, ".gz")
        open(filepath, "r") do file
            stream = GzipDecompressorStream(file)
            JSON3.read(stream)
        end
    else
        JSON3.read(read(filepath, String))
    end

    # Parse header
    header = parse_header(json_data.header)

    # Get all positions
    positions_json = json_data.positions
    total_positions = length(positions_json)

    if total_positions == 0
        @warn "No variants found in input file"
        return header
    end

    # Batch all variant indices for each thread
    chunk_size = ceil(Int, total_positions / nthreads)
    chunks = [collect(i:min(i+chunk_size-1, total_positions))
              for i in 1:chunk_size:total_positions]

    # Create progress bar
    progress = Progress(total_positions, desc="Processing variants: ",
                       barlen=50, showspeed=true)

    # For thread-safe progress updates
    progress_lock = ReentrantLock()

    # Process batches in parallel
    Threads.@threads for chunk_idx in 1:length(chunks)
        thread_id = Threads.threadid()
        chunk = chunks[chunk_idx]
        stats = thread_stats[thread_id]

        for idx in chunk
            pos = positions_json[idx]

            # Prefilter check (while tracking statistics)
            prefilter_result = check_prefilter_with_stats(pos, config)

            if !prefilter_result.passed
                # Update statistics
                if prefilter_result.reason == :failed_depth
                    stats[:failed_depth] += 1
                elseif prefilter_result.reason == :failed_vaf
                    stats[:failed_vaf] += 1
                elseif prefilter_result.reason == :failed_af
                    stats[:failed_af] += 1
                end

                # Thread-safe progress update
                lock(progress_lock) do
                    next!(progress)
                end
                continue
            end

            # Parse variant
            variant = parse_position(pos)
            if variant === nothing
                lock(progress_lock) do
                    next!(progress)
                end
                continue
            end

            # Call callback (pass thread ID)
            callback(variant, idx, thread_id)

            # Thread-safe progress update
            lock(progress_lock) do
                next!(progress)
            end
        end
    end

    # Finish progress bar
    finish!(progress)

    return header
end

"""
    PrefilterResult

Prefilter result
"""
struct PrefilterResult
    passed::Bool
    reason::Union{Symbol, Nothing}
end

"""
    check_prefilter_with_stats(pos_json, config) -> PrefilterResult

Quickly check if variant passes prefilter and return failure reason (for statistics tracking)

# Returns
- `PrefilterResult(true, nothing)`: Passed
- `PrefilterResult(false, :failed_depth)`: Insufficient depth
- `PrefilterResult(false, :failed_vaf)`: VAF too low
- `PrefilterResult(false, :failed_af)`: Population frequency too high
"""
function check_prefilter_with_stats(pos_json, config::FilterConfig)::PrefilterResult
    # 1. Check sequencing quality
    if haskey(pos_json, :samples) && !isempty(pos_json.samples)
        sample = pos_json.samples[1]

        # Check depth
        if haskey(sample, :totalDepth)
            depth = sample.totalDepth
            if depth < config.min_total_depth
                return PrefilterResult(false, :failed_depth)
            end
        end

        # Check VAF
        if haskey(sample, :variantFrequencies) && !isempty(sample.variantFrequencies)
            vaf = sample.variantFrequencies[1]
            if vaf < config.min_variant_frequency
                return PrefilterResult(false, :failed_vaf)
            end
        end
    end

    # 2. Check population frequency
    if haskey(pos_json, :variants) && !isempty(pos_json.variants)
        variant = pos_json.variants[1]

        # Check gnomad-exome EAS AF
        if haskey(variant, Symbol("gnomad-exome"))
            gne = variant[Symbol("gnomad-exome")]
            if haskey(gne, :easAf) && gne.easAf > config.max_eas_af
                return PrefilterResult(false, :failed_af)
            end
        end

        # Check 1000 Genomes EAS AF
        if haskey(variant, :oneKg)
            okg = variant.oneKg
            if haskey(okg, :easAf) && okg.easAf > config.max_eas_af
                return PrefilterResult(false, :failed_af)
            end
        end
    end

    # Pass all checks
    return PrefilterResult(true, nothing)
end

end # module NirvanaParser
