use anyhow::{Context, Result};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use json2maf::*;
use rayon::prelude::*;
use std::fs;
use std::path::Path;

#[derive(Parser, Debug)]
#[command(name = "json2maf")]
#[command(author = "JSON2MAF Contributors")]
#[command(version = "0.5.0")]
#[command(about = "Pathogenic variant filtering tool for Nirvana JSON", long_about = None)]
struct Args {
    /// Input Nirvana JSON.gz file path
    #[arg(short, long)]
    input: String,

    /// Output MAF file path
    #[arg(short, long)]
    output: String,

    /// Minimum sequencing depth
    #[arg(long, default_value_t = 30)]
    min_depth: i32,

    /// Minimum variant allele frequency VAF
    #[arg(long, default_value_t = 0.03)]
    min_vaf: f64,

    /// Maximum East Asian population allele frequency
    #[arg(long, default_value_t = 0.01)]
    max_eas_af: f64,

    /// REVEL score threshold
    #[arg(long, default_value_t = 0.75)]
    min_revel: f64,

    /// PrimateAI-3D score threshold
    #[arg(long, default_value_t = 0.8)]
    min_primate_ai: f64,

    /// DANN score threshold
    #[arg(long, default_value_t = 0.96)]
    min_dann: f64,

    /// Keep temporary files
    #[arg(long)]
    keep_temp: bool,

    /// Statistics report output path
    #[arg(long)]
    stats: Option<String>,

    /// Verbose output mode
    #[arg(short, long)]
    verbose: bool,

    /// Quiet mode (no progress display)
    #[arg(short, long)]
    quiet: bool,

    /// Number of threads (defaults to number of CPU cores)
    #[arg(short = 'j', long)]
    threads: Option<usize>,

    /// Exclude benign and likely benign variants
    #[arg(long)]
    exclude_benign: bool,
}

fn main() -> Result<()> {
    env_logger::Builder::from_default_env()
        .filter_level(log::LevelFilter::Info)
        .init();

    let args = Args::parse();

    // Set thread pool size
    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .context("Failed to set thread pool size")?;
    }

    let num_threads = rayon::current_num_threads();

    // Create configuration
    let config = FilterConfig {
        min_total_depth: args.min_depth,
        min_variant_frequency: args.min_vaf,
        max_eas_af: args.max_eas_af,
        min_revel_score: args.min_revel,
        min_primate_ai_score: args.min_primate_ai,
        min_dann_score: args.min_dann,
        exclude_benign: args.exclude_benign,
    };

    // Validate configuration
    config.validate()?;

    // Check input file exists
    if !Path::new(&args.input).exists() {
        anyhow::bail!("Input file does not exist: {}", args.input);
    }

    if args.verbose {
        println!("\nStarting processing: {}", args.input);
        println!("Output file: {}", args.output);
        println!("Processing mode: Multi-threaded parallel processing");
        println!("Number of threads: {}", num_threads);
        display_config(&config);
    }

    // Process file
    let stats = process_nirvana_json(
        &args.input,
        &args.output,
        &config,
        args.verbose,
        args.quiet,
        args.keep_temp,
    )?;

    // Print statistics
    if args.verbose || args.stats.is_some() {
        print_statistics(&stats, num_threads, args.stats.as_deref())?;
    }

    println!("\n✓ Processing complete! (Using {} threads for parallel processing)", num_threads);

    Ok(())
}

fn display_config(config: &FilterConfig) {
    println!("============================================================");
    println!("JSON2MAF Filter Configuration");
    println!("============================================================");
    println!();
    println!("Quality filtering parameters:");
    println!("  Minimum sequencing depth (min_total_depth):       {}", config.min_total_depth);
    println!("  Minimum VAF (min_variant_frequency):              {}", config.min_variant_frequency);
    println!();
    println!("Population frequency filtering parameters:");
    println!("  Maximum East Asian AF (max_eas_af):               {}", config.max_eas_af);
    println!();
    println!("Predictive score thresholds:");
    println!("  REVEL minimum score (min_revel_score):            {}", config.min_revel_score);
    println!("  PrimateAI-3D minimum score:                       {}", config.min_primate_ai_score);
    println!("  DANN minimum score:                               {}", config.min_dann_score);
    println!();
    println!("ClinVar filtering options:");
    println!("  Exclude benign/likely benign variants:            {}", config.exclude_benign);
    println!();
    println!("============================================================");
}

fn process_nirvana_json(
    input_path: &str,
    output_path: &str,
    config: &FilterConfig,
    verbose: bool,
    quiet: bool,
    _keep_temp: bool,
) -> Result<FilterStats> {
    if verbose {
        println!("\n[1/3] Parsing Nirvana JSON...");
    }

    // Parse JSON
    let (_header, variants) = parse_nirvana_json(input_path)?;

    if verbose {
        println!("  ✓ Parsed {} variants", variants.len());
        println!("\n[2/3] Filtering variants in parallel...");
    }

    let progress = if !quiet {
        let pb = ProgressBar::new(variants.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
                .unwrap()
                .progress_chars("#>-"),
        );
        Some(pb)
    } else {
        None
    };

    // Process variants in parallel
    let results: Vec<_> = variants
        .par_iter()
        .map(|variant| {
            let mut thread_stats = FilterStats::default();

            // Quality filtering
            let quality_result = apply_quality_filters(variant, config);

            if !quality_result.passes_quality {
                if let Some(reason) = &quality_result.failure_reason {
                    let reason_lower = reason.to_lowercase();
                    if reason_lower.contains("depth") {
                        thread_stats.failed_depth += 1;
                    } else if reason_lower.contains("frequency") || reason_lower.contains("vaf") {
                        thread_stats.failed_vaf += 1;
                    } else if reason_lower.contains("population")
                        || reason_lower.contains("allele frequency")
                        || reason_lower.contains("east asian")
                    {
                        thread_stats.failed_af += 1;
                    }
                }
                return (None, thread_stats);
            }

            thread_stats.passed_quality += 1;

            // ClinVar assessment
            let clinvar_assessment = assess_clinvar_pathogenicity(&variant.clinvar);

            // Predictive scores assessment
            let predictive_assessment = assess_predictive_scores(variant, config);

            // Integrated decision
            let decision = make_filter_decision_with_config(
                variant,
                &clinvar_assessment,
                &predictive_assessment,
                config.exclude_benign,
            );

            // Update statistics
            if decision.should_include {
                thread_stats.included += 1;

                if decision.primary_evidence == "ClinVar" {
                    if decision.pathogenicity_class == "Pathogenic" {
                        thread_stats.clinvar_pathogenic += 1;
                    } else if decision.pathogenicity_class == "Likely pathogenic" {
                        thread_stats.clinvar_likely += 1;
                    }
                } else if decision.primary_evidence == "Predictive" {
                    thread_stats.predictive_likely += 1;
                    if has_primate_ai_support(&predictive_assessment)
                        && count_supporting_predictive_scores(&predictive_assessment) == 1
                    {
                        thread_stats.primate_ai_only += 1;
                    } else {
                        thread_stats.multi_score += 1;
                    }
                }

                // Convert to MAF record
                let maf_record = variant_to_maf(variant, &decision);

                if let Some(pb) = &progress {
                    pb.inc(1);
                }

                (Some(maf_record), thread_stats)
            } else {
                thread_stats.excluded += 1;

                // Track benign exclusions separately
                if decision.pathogenicity_class.contains("Benign") {
                    thread_stats.excluded_benign += 1;
                }

                if let Some(pb) = &progress {
                    pb.inc(1);
                }

                (None, thread_stats)
            }
        })
        .collect();

    if let Some(pb) = progress {
        pb.finish_with_message("Filtering complete");
    }

    // Merge results
    let mut total_stats = FilterStats::default();
    let mut maf_records = Vec::new();

    for (record, stats) in results {
        total_stats.merge(&stats);
        if let Some(rec) = record {
            maf_records.push(rec);
        }
    }

    if verbose {
        println!("  ✓ Filtered {} / {} variants", total_stats.included, variants.len());
        println!("\n[3/3] Writing MAF file...");
    }

    // Write MAF file
    let mut writer = MAFWriter::new(output_path)?;
    for record in &maf_records {
        writer.write_record(record)?;
    }
    writer.flush()?;

    if verbose {
        println!("  ✓ Successfully wrote {} records to {}", total_stats.included, output_path);
    }

    Ok(total_stats)
}

fn print_statistics(stats: &FilterStats, num_threads: usize, output_path: Option<&str>) -> Result<()> {
    let benign_section = if stats.excluded_benign > 0 {
        format!("\nClinVar benign filtering:\n  - Excluded benign/likely benign: {}\n", stats.excluded_benign)
    } else {
        String::new()
    };

    let report = format!(
        r#"
═══════════════════════════════════════════════════════════
                  Filtering Statistics Report
═══════════════════════════════════════════════════════════

Processing mode:        Multi-threaded parallel processing
Number of threads:      {}

Quality filtering:
  - Passed quality:     {}
  - Insufficient depth: {}
  - VAF too low:        {}
  - Population freq too high: {}

Pathogenicity assessment:
  - ClinVar Pathogenic:         {}
  - ClinVar Likely pathogenic:  {}
  - Predictive scores support:  {}
    * PrimateAI-3D solo support: {}
    * 2+ scores support:         {}
{}
Final results:
  - Included variants:  {}
  - Excluded variants:  {}

═══════════════════════════════════════════════════════════
"#,
        num_threads,
        stats.passed_quality,
        stats.failed_depth,
        stats.failed_vaf,
        stats.failed_af,
        stats.clinvar_pathogenic,
        stats.clinvar_likely,
        stats.predictive_likely,
        stats.primate_ai_only,
        stats.multi_score,
        benign_section,
        stats.included,
        stats.excluded
    );

    println!("{}", report);

    if let Some(path) = output_path {
        fs::write(path, report).context("Failed to write statistics report")?;
        println!("Statistics report written to: {}", path);
    }

    Ok(())
}
