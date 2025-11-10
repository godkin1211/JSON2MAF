# JSON2MAF - Rust Version

A high-performance Rust implementation of JSON2MAF for filtering pathogenic variants from Illumina Nirvana annotation output and converting them to MAF (Mutation Annotation Format).

[![Rust](https://img.shields.io/badge/Rust-1.70+-orange.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/License-Internal-red.svg)]()

## Overview

This Rust version provides the same functionality as the Julia implementation with enhanced performance and memory efficiency. It processes variant annotation data from Illumina Nirvana (gzipped JSON format) and extracts pathogenic/likely pathogenic variants for oncology research.

### Key Features

- **Multi-threaded parallel processing** - Utilizes Rayon for efficient parallel processing
- **Quality filtering** - Filters variants based on sequencing depth and variant allele frequency (VAF)
- **Population frequency filtering** - Excludes common variants in East Asian populations
- **ClinVar prioritization** - Intelligent handling of multiple ClinVar entries with conflict resolution
- **Predictive score integration** - Combines REVEL, DANN, PrimateAI-3D, and COSMIC scores
- **Complete statistics tracking** - Provides detailed reports on filtering stages
- **Configurable thresholds** - All filtering parameters can be customized
- **Zero-cost abstractions** - Rust's performance guarantees for production use

### Performance Benefits

The Rust implementation offers several advantages:

- **Memory safety** - No runtime overhead with compile-time guarantees
- **Zero-copy parsing** - Efficient JSON processing with serde
- **Optimized binary** - Aggressive compilation optimizations enabled
- **Lower memory footprint** - Efficient data structures and memory management

## System Requirements

- **Rust**: 1.70 or higher
- **RAM**: 2-4 GB (varies with input file size and thread count)
- **CPU**: Multi-core processor recommended for optimal performance

## Installation

### Building from Source

```bash
# Navigate to the Rust version directory
cd rust_version

# Build in release mode (optimized)
cargo build --release

# The binary will be available at target/release/json2maf
```

### Running Tests

```bash
cargo test
```

## Usage

### Quick Start

```bash
# Basic usage
./target/release/json2maf \
  -i input.json.gz \
  -o output.maf

# With verbose output and statistics report
./target/release/json2maf \
  -i input.json.gz \
  -o output.maf \
  --stats report.txt \
  --verbose

# Specify number of threads
./target/release/json2maf \
  -i input.json.gz \
  -o output.maf \
  -j 8
```

### Complete Example with Custom Parameters

```bash
./target/release/json2maf \
  --input input.json.gz \
  --output output.maf \
  --min-depth 30 \
  --min-vaf 0.03 \
  --max-eas-af 0.01 \
  --min-revel 0.75 \
  --min-primate-ai 0.8 \
  --min-dann 0.96 \
  --threads 8 \
  --stats report.txt \
  --verbose
```

### Command-Line Options

**Required Arguments**:

- `-i, --input <FILE>`: Input Nirvana JSON.gz file path
- `-o, --output <FILE>`: Output MAF file path

**Quality Filtering Parameters**:

- `--min-depth <INT>`: Minimum sequencing depth (default: 30)
- `--min-vaf <FLOAT>`: Minimum variant allele frequency (default: 0.03)
- `--max-eas-af <FLOAT>`: Maximum East Asian allele frequency (default: 0.01)

**Predictive Score Thresholds**:

- `--min-revel <FLOAT>`: REVEL score threshold (default: 0.75)
- `--min-primate-ai <FLOAT>`: PrimateAI-3D score threshold (default: 0.8)
- `--min-dann <FLOAT>`: DANN score threshold (default: 0.96)

**Other Options**:

- `-j, --threads <NUM>`: Number of threads (default: number of CPU cores)
- `-v, --verbose`: Enable verbose output
- `-q, --quiet`: Suppress progress bar
- `--stats <FILE>`: Save statistics report to file
- `--keep-temp`: Keep temporary thread files (for debugging)
- `-h, --help`: Show help message
- `--version`: Show version information

### Performance Tuning

The Rust version automatically uses all available CPU cores by default. You can adjust this:

```bash
# Use 8 threads
./target/release/json2maf -i input.json.gz -o output.maf -j 8

# Use all available cores (default)
./target/release/json2maf -i input.json.gz -o output.maf
```

**Recommended thread counts**:

- 4 threads: Standard workstations
- 6-8 threads: High-performance workstations
- 16+ threads: Server environments

## Output Format

### MAF File

The tool generates a standard MAF (Mutation Annotation Format) file with the same fields as the Julia version. See the main README for detailed field descriptions.

### Statistics Report

```
═══════════════════════════════════════════════════════════
                  Filtering Statistics Report
═══════════════════════════════════════════════════════════

Processing mode:        Multi-threaded parallel processing
Number of threads:      8

Quality filtering:
  - Passed quality:     1,340
  - Insufficient depth: 623
  - VAF too low:        1,193
  - Population freq too high: 27,545

Pathogenicity assessment:
  - ClinVar Pathogenic:         0
  - ClinVar Likely pathogenic:  0
  - Predictive scores support:  198
    * PrimateAI-3D solo support: 7
    * 2+ scores support:         191

Final results:
  - Included variants:  198
  - Excluded variants:  1,142

═══════════════════════════════════════════════════════════
```

## Project Structure

```
rust_version/
├── Cargo.toml              # Project configuration
├── src/
│   ├── main.rs             # Main executable with CLI
│   ├── lib.rs              # Library exports
│   ├── types.rs            # Data structures
│   ├── parser.rs           # JSON parsing
│   ├── filters/
│   │   ├── mod.rs          # Filter module exports
│   │   ├── quality.rs      # Quality filtering
│   │   ├── clinvar.rs      # ClinVar filtering
│   │   ├── predictive.rs   # Predictive scores
│   │   └── decision.rs     # Decision logic
│   ├── converter.rs        # MAF conversion
│   └── writer.rs           # MAF output
└── README.md               # This file
```

## Development

### Running Tests

```bash
cargo test
```

### Running with Logging

```bash
RUST_LOG=info ./target/release/json2maf -i input.json.gz -o output.maf
```

### Building for Maximum Performance

The release profile is already optimized for maximum performance:

```toml
[profile.release]
opt-level = 3         # Maximum optimization
lto = true            # Link-time optimization
codegen-units = 1     # Better optimization at cost of compile time
strip = true          # Strip symbols from binary
```

### Benchmarking

```bash
cargo bench
```

## Differences from Julia Version

1. **Dependency Management**: Uses Cargo instead of Julia's Pkg
2. **Parallel Processing**: Uses Rayon instead of Julia's native threading
3. **Error Handling**: Uses Result types instead of exceptions
4. **Type System**: Strongly typed with compile-time guarantees
5. **Performance**: Generally comparable or faster due to LLVM optimizations

## Troubleshooting

### Common Issues

**Issue**: Compilation errors

- **Solution**: Ensure Rust 1.70+ is installed: `rustc --version`

**Issue**: Out of memory errors

- **Solution**: Reduce thread count with `-j` flag

**Issue**: Slow compilation

- **Solution**: Use `cargo build --release` only for production builds

## Contributing

This is an internal research tool. For questions or suggestions, please contact:

- Project bioinformaticians
- Product managers

## License

Internal use only for oncology research and drug development.

## Acknowledgments

- Original Julia implementation by the JSON2MAF team
- Rust community for excellent libraries (serde, rayon, clap, etc.)

---

**Version**: 0.4.0
**Last Updated**: November 2025
**Status**: Production-ready
