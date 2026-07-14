use anyhow::{Context, Result};
use clap::Parser;
use json2maf::sv::{parse_sv_nirvana_streaming, sv_position_to_record, SVType, SVWriter};
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

    let mut writer = SVWriter::new(&args.output).context("Failed to create output file")?;

    let mut total = 0usize;
    let mut del_count = 0usize;
    let mut ins_count = 0usize;
    let mut symbolic = 0usize;

    let _header = parse_sv_nirvana_streaming(&args.input, |pos| {
        total += 1;
        match pos.sv_type {
            SVType::Del => del_count += 1,
            SVType::Ins => ins_count += 1,
        }
        if pos.is_symbolic {
            symbolic += 1;
        }

        let record = sv_position_to_record(&pos);
        writer.write_record(&record)
    })
    .context("Failed to parse SV JSON")?;

    writer.flush()?;

    if args.verbose {
        println!(
            "Parsed {} SV positions: {} DEL, {} INS ({} symbolic, {} sequence-resolved)",
            total,
            del_count,
            ins_count,
            symbolic,
            total - symbolic,
        );
    }

    println!("Written {} SV records to {}", total, args.output);

    Ok(())
}
