use anyhow::{Context, Result};
use clap::Parser;
use json2maf::sv::{parse_sv_nirvana_json, sv_position_to_record, SVWriter};
use std::path::Path;

#[derive(Parser, Debug)]
#[command(name = "json2sv")]
#[command(version = "0.1.0")]
#[command(about = "Convert Nirvana SV-annotated JSON to TSV (INS/DEL only)", long_about = None)]
struct Args {
    /// Input Nirvana SV JSON.gz file
    #[arg(short, long)]
    input: String,

    /// Output TSV file path
    #[arg(short, long)]
    output: String,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

fn main() -> Result<()> {
    env_logger::Builder::from_default_env()
        .filter_level(log::LevelFilter::Warn)
        .init();

    let args = Args::parse();

    if !Path::new(&args.input).exists() {
        anyhow::bail!("Input file not found: {}", args.input);
    }

    if args.verbose {
        println!("Input:  {}", args.input);
        println!("Output: {}", args.output);
    }

    let (_header, sv_positions) =
        parse_sv_nirvana_json(&args.input).context("Failed to parse SV JSON")?;

    if args.verbose {
        let del_count = sv_positions.iter().filter(|p| p.sv_type == json2maf::sv::SVType::Del).count();
        let ins_count = sv_positions.iter().filter(|p| p.sv_type == json2maf::sv::SVType::Ins).count();
        let symbolic = sv_positions.iter().filter(|p| p.is_symbolic).count();
        println!(
            "Parsed {} SV positions: {} DEL, {} INS ({} symbolic, {} sequence-resolved)",
            sv_positions.len(),
            del_count,
            ins_count,
            symbolic,
            sv_positions.len() - symbolic,
        );
    }

    let mut writer = SVWriter::new(&args.output).context("Failed to create output file")?;

    for pos in &sv_positions {
        let record = sv_position_to_record(pos);
        writer.write_record(&record)?;
    }
    writer.flush()?;

    println!(
        "Written {} SV records to {}",
        sv_positions.len(),
        args.output
    );

    Ok(())
}
