#!/usr/bin/env julia

"""
json2maf.jl

Pathogenic variant filtering tool for Nirvana JSON

Multi-threaded parallel processing for large-scale Nirvana JSON files:
- Utilize multi-core CPU for parallel variant processing
- Complete statistics tracking (including pre-filtering stage)
- Each thread processes a batch of variants independently and merges results
- Supports high-performance processing of large VCF annotation files

⚠️ Important: Recommended to set thread count for optimal performance
export JULIA_NUM_THREADS=8

Usage:
    julia -t 8 bin/json2maf.jl -i input.json.gz -o output.maf [options]
"""

using Pkg
# Activate project environment
Pkg.activate(joinpath(@__DIR__, ".."))

using ArgParse
using JSON2MAF

"""
    parse_commandline()

Parse command line arguments
"""
function parse_commandline()
    s = ArgParseSettings(
        prog = "json2maf",
        description = "Pathogenic variant filtering tool for Nirvana JSON",
        version = string(JSON2MAF.version()),
        add_version = true
    )

    @add_arg_table! s begin
        "--input", "-i"
            help = "Input Nirvana JSON.gz file path"
            required = true
            arg_type = String

        "--output", "-o"
            help = "Output MAF file path"
            required = true
            arg_type = String

        "--min-depth"
            help = "Minimum sequencing depth (default: 30)"
            arg_type = Int
            default = 30

        "--min-vaf"
            help = "Minimum variant allele frequency VAF (default: 0.03)"
            arg_type = Float64
            default = 0.03

        "--max-eas-af"
            help = "Maximum East Asian population allele frequency (default: 0.01)"
            arg_type = Float64
            default = 0.01

        "--min-revel"
            help = "REVEL score threshold (default: 0.75)"
            arg_type = Float64
            default = 0.75

        "--min-primate-ai"
            help = "PrimateAI-3D score threshold (default: 0.8)"
            arg_type = Float64
            default = 0.8

        "--min-dann"
            help = "DANN score threshold (default: 0.96)"
            arg_type = Float64
            default = 0.96

        "--keep-temp"
            help = "Keep temporary files (default: delete)"
            action = :store_true

        "--stats"
            help = "Statistics report output path (optional)"
            arg_type = String
            default = nothing

        "--verbose", "-v"
            help = "Verbose output mode"
            action = :store_true

        "--quiet", "-q"
            help = "Quiet mode (no progress display)"
            action = :store_true
    end

    return parse_args(s)
end

"""
    print_statistics(stats::Dict, output_path::Union{String, Nothing}=nothing)

Output statistics
"""
function print_statistics(stats::Dict, output_path::Union{String, Nothing}=nothing)
    report = """

    ═══════════════════════════════════════════════════════════
                      Filtering Statistics Report
    ═══════════════════════════════════════════════════════════

    Processing mode:        Multi-threaded parallel processing
    Number of threads:      $(stats[:num_threads])

    Quality filtering:
      - Passed quality:     $(stats[:passed_quality])
      - Insufficient depth: $(stats[:failed_depth])
      - VAF too low:        $(stats[:failed_vaf])
      - Population freq too high: $(stats[:failed_af])

    Pathogenicity assessment:
      - ClinVar Pathogenic:         $(stats[:clinvar_pathogenic])
      - ClinVar Likely pathogenic:  $(stats[:clinvar_likely])
      - Predictive scores support:  $(stats[:predictive_likely])
        * PrimateAI-3D solo support: $(stats[:primate_ai_only])
        * 2+ scores support:         $(stats[:multi_score])

    Final results:
      - Included variants:  $(stats[:included])
      - Excluded variants:  $(stats[:excluded])

    ═══════════════════════════════════════════════════════════
    """

    println(report)

    # If statistics file path is specified, write to file
    if output_path !== nothing
        open(output_path, "w") do io
            write(io, report)
        end
        println("Statistics report written to: $output_path")
    end
end

"""
    process_nirvana_json(input_path::String, output_path::String, config::FilterConfig;
                        keep_temp::Bool=false, verbose::Bool=false, quiet::Bool=false,
                        stats_path::Union{String,Nothing}=nothing)

Process Nirvana JSON file and output as MAF format
"""
function process_nirvana_json(input_path::String, output_path::String, config::FilterConfig;
                                      keep_temp::Bool=false, verbose::Bool=false, quiet::Bool=false,
                                      stats_path::Union{String,Nothing}=nothing)

    nthreads = Threads.nthreads()

    if verbose
        println("\nStarting processing: $input_path")
        println("Output file: $output_path")
        println("Processing mode: Multi-threaded parallel processing")
        println("Number of threads: $nthreads")
        display_config(config)
    end

    if nthreads == 1
        @warn "Only 1 thread! Recommended to set JULIA_NUM_THREADS for better performance."
        @warn "Example: export JULIA_NUM_THREADS=8"
    end

    # Create independent statistics for each thread
    thread_stats = [Dict(
        :passed_quality => 0,
        :failed_depth => 0,
        :failed_vaf => 0,
        :failed_af => 0,
        :clinvar_pathogenic => 0,
        :clinvar_likely => 0,
        :predictive_likely => 0,
        :primate_ai_only => 0,
        :multi_score => 0,
        :included => 0,
        :excluded => 0
    ) for _ in 1:nthreads]

    # Create independent MAF writer for each thread
    temp_files = ["$(output_path).thread_$(i).tmp" for i in 1:nthreads]
    thread_writers = [create_maf_writer(temp_files[i], batch_size=1000) for i in 1:nthreads]

    if verbose
        println("\n[1/4] Parallel processing Nirvana JSON...")
    end

    # Process all variants in parallel (track complete statistics from pre-filtering stage)
    try
        header = process_nirvana_parallel_with_stats(input_path, config, thread_stats) do variant, idx, tid
            stats = thread_stats[tid]
            writer = thread_writers[tid]

            # Quality filtering
            quality_result = apply_quality_filters(variant, config)

            if !quality_result.passes_quality
                # Update statistics - categorize by failure reason
                if quality_result.failure_reason !== nothing
                    reason_lower = lowercase(quality_result.failure_reason)
                    if occursin("depth", reason_lower)
                        stats[:failed_depth] += 1
                    elseif occursin("frequency", reason_lower) || occursin("vaf", reason_lower)
                        stats[:failed_vaf] += 1
                    elseif occursin("population", reason_lower) ||
                           occursin("allele frequency", reason_lower) ||
                           occursin("east asian", reason_lower)
                        stats[:failed_af] += 1
                    end
                end
                return
            end

            stats[:passed_quality] += 1

            # ClinVar assessment
            clinvar_assessment = assess_clinvar_pathogenicity(variant.clinvar)

            # Predictive scores assessment
            predictive_assessment = assess_predictive_scores(variant, config)

            # Integrated decision
            decision = make_filter_decision(variant, clinvar_assessment, predictive_assessment)

            # Update statistics and write
            if decision.should_include
                stats[:included] += 1

                if decision.primary_evidence == "ClinVar"
                    if decision.pathogenicity_class == "Pathogenic"
                        stats[:clinvar_pathogenic] += 1
                    elseif decision.pathogenicity_class == "Likely pathogenic"
                        stats[:clinvar_likely] += 1
                    end
                elseif decision.primary_evidence == "Predictive"
                    stats[:predictive_likely] += 1
                    if has_primate_ai_support(predictive_assessment) &&
                       count_supporting_predictive_scores(predictive_assessment) == 1
                        stats[:primate_ai_only] += 1
                    else
                        stats[:multi_score] += 1
                    end
                end

                # Convert to MAF record and write
                maf_record = variant_to_maf(variant, decision)
                write_maf_batch(writer, maf_record)

                # Batch output (output every 50 variants to reduce I/O overhead)
                if verbose && stats[:included] % 50 == 0
                    println("  ✓ [Thread $tid] Processed $(stats[:included]) variants " *
                           "(latest: $(maf_record.hugo_symbol) $(maf_record.chromosome):$(maf_record.start_position))")
                end
            else
                stats[:excluded] += 1
            end
        end

    finally
        # Close all writers
        for writer in thread_writers
            close_maf_writer(writer)
        end
    end

    if verbose
        println("\n[2/4] Merging thread outputs...")
    end

    # Merge all thread output files
    merge_maf_files(temp_files, output_path, keep_temp=keep_temp)

    # Merge statistics
    total_stats = Dict(
        :num_threads => nthreads,
        :passed_quality => sum(s[:passed_quality] for s in thread_stats),
        :failed_depth => sum(s[:failed_depth] for s in thread_stats),
        :failed_vaf => sum(s[:failed_vaf] for s in thread_stats),
        :failed_af => sum(s[:failed_af] for s in thread_stats),
        :clinvar_pathogenic => sum(s[:clinvar_pathogenic] for s in thread_stats),
        :clinvar_likely => sum(s[:clinvar_likely] for s in thread_stats),
        :predictive_likely => sum(s[:predictive_likely] for s in thread_stats),
        :primate_ai_only => sum(s[:primate_ai_only] for s in thread_stats),
        :multi_score => sum(s[:multi_score] for s in thread_stats),
        :included => sum(s[:included] for s in thread_stats),
        :excluded => sum(s[:excluded] for s in thread_stats)
    )

    total_written = count_maf_records(output_path)

    if verbose
        println("\n[3/4] Completed MAF file writing...")
        println("  ✓ Successfully merged and wrote $total_written records to $output_path")
        if !keep_temp
            println("  ✓ Deleted temporary files")
        end
    end

    # Output statistics
    if verbose
        println("\n[4/4] Generating statistics report...")
    end

    print_statistics(total_stats, stats_path)

    return total_stats
end

"""
    main()

Main program entry point
"""
function main()
    # Parse arguments
    args = parse_commandline()

    # Create configuration
    config = create_filter_config(
        min_total_depth = args["min-depth"],
        min_variant_frequency = args["min-vaf"],
        max_eas_af = args["max-eas-af"],
        min_revel_score = args["min-revel"],
        min_primate_ai_score = args["min-primate-ai"],
        min_dann_score = args["min-dann"]
    )

    # Validate configuration
    try
        validate_config(config)
    catch e
        println(stderr, "Error: Invalid configuration parameters - $e")
        exit(1)
    end

    # Check input file
    if !isfile(args["input"])
        println(stderr, "Error: Input file does not exist: $(args["input"])")
        exit(1)
    end

    # Check thread count
    nthreads = Threads.nthreads()
    if nthreads == 1 && !args["quiet"]
        println("Warning: Only 1 thread!")
        println("         Recommended to set JULIA_NUM_THREADS environment variable to enable multi-threading")
        println("         Example: export JULIA_NUM_THREADS=8")
        println("         Or use: julia -t 8 $(PROGRAM_FILE)")
        println()
    end

    try
        # Process file
        process_nirvana_json(
            args["input"],
            args["output"],
            config,
            keep_temp = args["keep-temp"],
            verbose = args["verbose"],
            quiet = args["quiet"],
            stats_path = args["stats"]
        )

        println("\n✓ Processing complete! (Using $nthreads threads for parallel processing)")
        exit(0)

    catch e
        println(stderr, "\n✗ Error occurred: $e")
        if args["verbose"]
            showerror(stderr, e, catch_backtrace())
        end
        exit(1)
    end
end

# Execute main program
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
