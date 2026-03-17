# JSON2MAF - Rust Version

A high-performance Rust implementation for filtering pathogenic variants from Illumina Nirvana annotation output and converting them to MAF (Mutation Annotation Format).

[![Rust](https://img.shields.io/badge/Rust-1.70+-orange.svg)](https://www.rust-lang.org/)
[![Tests](https://img.shields.io/badge/Tests-25%20passing-brightgreen.svg)]()
[![License](https://img.shields.io/badge/License-Internal-red.svg)]()

## Overview

This Rust implementation provides identical functionality to the Julia version with enhanced performance, memory safety, and type guarantees. It processes variant annotation data from Illumina Nirvana (gzipped JSON format) and extracts pathogenic/likely pathogenic variants for oncology research.

### Key Features

- **Complete MAF annotation** - 37 output fields including transcript details, consequence predictions, and clinical significance
- **Multi-threaded parallel processing** - Utilizes Rayon for efficient parallel variant processing
- **Intelligent filtering pipeline** - Multi-stage quality, population frequency, and pathogenicity assessment
- **ClinVar prioritization** - Conflict resolution with cancer-specific prioritization
- **Predictive score integration** - REVEL, DANN, PrimateAI-3D with evidence-based thresholds
- **Comprehensive statistics** - Detailed filtering reports for quality control
- **Memory safety** - Compile-time guarantees with zero runtime overhead
- **Configurable thresholds** - All filtering parameters can be customized via CLI

### Performance Benefits

- **Zero-copy parsing** - Efficient JSON processing with serde
- **Parallel processing** - Rayon-based multi-threading with automatic load balancing
- **Optimized binary** - LTO and aggressive optimizations enabled
- **Unique temporary files** - Thread-safe with automatic cleanup
- **Lower memory footprint** - Efficient data structures (~2-4 GB for typical workloads)

## System Requirements

- **Rust**: 1.70 or higher
- **RAM**: 2-4 GB (scales with input size and thread count)
- **CPU**: Multi-core processor recommended (4-8 cores optimal)
- **Disk**: Temporary space for decompression (~2x input file size)

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
# Run all tests (25 tests)
cargo test

# Run with output
cargo test -- --nocapture

# Run specific test
cargo test test_parse_and_convert_with_new_fields
```

**Test Coverage**: 25 tests covering:
- JSON parsing and field extraction
- Quality filtering logic
- ClinVar assessment and conflict resolution
- Predictive score evaluation
- Decision engine logic
- MAF conversion and formatting
- Integration tests with real-world data

## Usage

### Quick Start

```bash
# Basic usage
./target/release/json2maf \
  -i input.json.gz \
  -o output.maf

# With statistics and verbose output
./target/release/json2maf \
  -i input.json.gz \
  -o output.maf \
  --stats report.txt \
  --verbose

# Specify number of threads
./target/release/json2maf \
  -i input.json.gz \
  -o output.maf \
  -j 6
```

### Complete Example with Custom Parameters

```bash
./target/release/json2maf \
  --input P.hard-filtered.vcf.annotated.json.gz \
  --output pathogenic_variants.maf \
  --min-depth 30 \
  --min-vaf 0.03 \
  --max-eas-af 0.01 \
  --min-revel 0.75 \
  --min-primate-ai 0.8 \
  --min-dann 0.96 \
  --threads 6 \
  --stats filtering_report.txt \
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

**Performance Options**:

- `-j, --threads <NUM>`: Number of threads (default: number of CPU cores)
- `-v, --verbose`: Enable verbose output
- `-q, --quiet`: Suppress progress bar
- `--keep-temp`: Keep temporary thread files (for debugging)

**Output Options**:

- `--stats <FILE>`: Save detailed statistics report to file

**Other**:

- `-h, --help`: Show help message
- `--version`: Show version information

## Filtering Logic

The tool implements a hierarchical decision engine:

1. **Quality Filtering** → Depth ≥30, VAF ≥0.03, EAS AF ≤0.01
2. **ClinVar Pathogenic** → Include as "Pathogenic"
3. **ClinVar Likely Pathogenic** → Include as "Likely pathogenic"
4. **ClinVar Inconclusive + Predictive Support** → Include as "Likely pathogenic"
   - PrimateAI-3D alone (threshold 0.8), OR
   - 2+ scores from {REVEL ≥0.75, DANN ≥0.96, PrimateAI-3D ≥0.8}
5. **All Other Cases** → Exclude

### ClinVar Conflict Resolution

When multiple ClinVar entries conflict:
- Prioritize cancer-related phenotypes
- Higher review status takes precedence
- Pathogenic > Likely pathogenic > Inconclusive

## Output Format

### MAF File (37 columns)

The tool generates a standard MAF file with the following fields:

**Basic Information**:
- `Hugo_Symbol` - Gene symbol (HGNC)
- `Chromosome` - Chromosome name
- `Start_Position` - Variant start position
- `End_Position` - Variant end position
- `Strand` - Genomic strand (always "+")

**Variant Classification**:
- `Variant_Classification` - MAF classification (Missense_Mutation, Nonsense_Mutation, etc.)
- `Variant_Type` - Type (SNP, INS, DEL, DNP, etc.)
- `Reference_Allele` - Reference allele
- `Tumor_Seq_Allele1` - Tumor allele 1 (usually same as reference)
- `Tumor_Seq_Allele2` - Tumor allele 2 (variant allele)

**Sample Information**:
- `Tumor_Sample_Barcode` - Sample identifier

**Transcript Annotation**:
- `HGVSc` - HGVS coding notation (e.g., "c.1799T>A")
- `HGVSp` - HGVS protein notation (e.g., "p.Val600Glu")
- `HGVSp_Short` - Short HGVS protein (e.g., "p.V600E")
- `Transcript_ID` - RefSeq or Ensembl transcript ID

**NEW: Detailed Consequence Annotation**:
- `Exon` - Exon number and total (e.g., "15/18")
- `Consequence` - Sequence Ontology terms (comma-separated)
- `IMPACT` - VEP-style impact level (HIGH, MODERATE, LOW, MODIFIER)
- `Codons` - Reference/variant codons (e.g., "Gtg/Gag")
- `Amino_Acids` - Reference/variant amino acids (e.g., "V/E")
- `cDNA_position` - Position in cDNA (e.g., "1799/2301")
- `CDS_position` - Position in CDS (e.g., "1799/2301")
- `Protein_position` - Position in protein (e.g., "600/766")

**Database IDs**:
- `dbSNP_RS` - dbSNP rs identifier
- `dbSNP_Val_Status` - Validation status (usually empty)
- `COSMIC_ID` - COSMIC mutation identifier

**ClinVar Information**:
- `ClinVar_ID` - ClinVar RCV identifier
- `ClinVar_Review_Status` - Review status (e.g., "criteria provided, multiple submitters")
- `ClinVar_Significance` - Clinical significance (Pathogenic, Likely pathogenic, etc.)
- `ClinVar_Disease` - Associated diseases/phenotypes

**Predictive Scores**:
- `PrimateAI_Score` - PrimateAI-3D pathogenicity score (0-1)
- `DANN_Score` - DANN pathogenicity score (0-1)
- `REVEL_Score` - REVEL pathogenicity score (0-1)

**Population Frequencies**:
- `gnomAD_AF` - gnomAD overall allele frequency
- `gnomAD_EAS_AF` - gnomAD East Asian allele frequency

**Sequencing Quality**:
- `Depth` - Total sequencing depth
- `VAF` - Variant allele frequency

### Statistics Report

```
═══════════════════════════════════════════════════════════
                  Filtering Statistics Report
═══════════════════════════════════════════════════════════

Processing mode:        Multi-threaded parallel processing
Number of threads:      2

Quality filtering:
  - Passed quality:     1,574
  - Insufficient depth: 1,503
  - VAF too low:        519
  - Population freq too high: 23,451

Pathogenicity assessment:
  - ClinVar Pathogenic:         2
  - ClinVar Likely pathogenic:  2
  - Predictive scores support:  213
    * PrimateAI-3D solo support: 7
    * 2+ scores support:         206

Final results:
  - Included variants:  217
  - Excluded variants:  1,357

═══════════════════════════════════════════════════════════
```

## Project Structure

```
rust_version/
├── Cargo.toml              # Project configuration and dependencies
├── README.md               # This file
├── src/
│   ├── main.rs             # CLI executable entry point
│   ├── lib.rs              # Library exports
│   ├── types.rs            # Core data structures (FilterConfig, VariantPosition, MAFRecord, etc.)
│   ├── parser.rs           # Nirvana JSON parsing with gzip decompression
│   ├── filters/
│   │   ├── mod.rs          # Filter module exports
│   │   ├── quality.rs      # Quality and population frequency filtering
│   │   ├── clinvar.rs      # ClinVar assessment and conflict resolution
│   │   ├── predictive.rs   # Predictive score evaluation (REVEL, DANN, PrimateAI-3D)
│   │   └── decision.rs     # Hierarchical decision engine
│   ├── converter.rs        # MAF format conversion
│   └── writer.rs           # Multi-threaded MAF file writing
└── tests/
    └── integration_test.rs # Integration tests for end-to-end validation
```

## Development

### Running Tests

```bash
# Run all tests
cargo test

# Run with verbose output
cargo test -- --nocapture

# Run specific test module
cargo test converter

# Run integration tests only
cargo test --test integration_test
```

### Code Coverage

- **Parser tests**: JSON parsing, field extraction, edge cases
- **Filter tests**: Quality, ClinVar, predictive scores, decision logic
- **Converter tests**: Variant classification mapping, HGVS notation, field extraction
- **Integration tests**: End-to-end processing with real-world data

### Running with Logging

```bash
# Info level
RUST_LOG=info ./target/release/json2maf -i input.json.gz -o output.maf

# Debug level
RUST_LOG=debug ./target/release/json2maf -i input.json.gz -o output.maf
```

### Building for Maximum Performance

The release profile is optimized for production use:

```toml
[profile.release]
opt-level = 3         # Maximum optimization level
lto = true            # Link-time optimization
codegen-units = 1     # Better optimization (slower compile)
strip = true          # Strip debug symbols from binary
```

Build command:
```bash
cargo build --release
```

### Performance Tuning

**Recommended Thread Counts**:
- **4 threads**: Standard workstations, small files (<100 MB)
- **6-8 threads**: High-performance workstations, medium files (100-500 MB)
- **16+ threads**: Server environments, large files (>500 MB)

**Memory Considerations**:
- Each thread creates a temporary MAF file
- Peak memory usage: ~2-4 GB for typical workloads
- Reduce threads if encountering memory issues

## Differences from Julia Version

| Aspect | Julia | Rust |
|--------|-------|------|
| **Dependency Management** | Pkg | Cargo |
| **Parallel Processing** | `Threads.@threads` | Rayon |
| **Error Handling** | Exceptions | `Result<T, E>` types |
| **Type System** | Dynamic with optional typing | Strong static typing |
| **Memory Safety** | GC-based | Ownership system (compile-time) |
| **Performance** | Fast JIT compilation | Ahead-of-time compilation |
| **Binary Size** | Requires Julia runtime | Standalone executable |
| **Test Framework** | `@testset` macros | Built-in `#[test]` |

## Troubleshooting

### Common Issues

**Issue**: Compilation errors

**Solution**:
```bash
# Ensure Rust 1.70+ is installed
rustc --version

# Update Rust if needed
rustup update stable
```

**Issue**: Out of memory errors

**Solution**:
```bash
# Reduce thread count
./target/release/json2maf -i input.json.gz -o output.maf -j 2
```

**Issue**: Slow first-time compilation

**Solution**: This is expected. Release builds with LTO take longer but produce optimized binaries. Use `cargo build` (without `--release`) for faster development builds.

**Issue**: Tests fail due to temporary file conflicts

**Solution**: Tests now use unique temporary files automatically. If issues persist, run tests sequentially:
```bash
cargo test -- --test-threads=1
```

**Issue**: "Invalid gzip header" errors

**Solution**: Ensure input file is properly gzipped:
```bash
# Check file type
file input.json.gz

# Re-compress if needed
gunzip input.json.gz
gzip input.json
```

## Example Workflow

```bash
# 1. Build the tool
cd rust_version
cargo build --release

# 2. Run on sample data
./target/release/json2maf \
  -i ../examples/sample_data/test_mini.json.gz \
  -o sample_output.maf \
  --verbose

# 3. Run on real data with statistics
./target/release/json2maf \
  -i ../P.hard-filtered.vcf.annotated.json.gz \
  -o pathogenic_variants.maf \
  --stats filtering_report.txt \
  -j 6 \
  --verbose

# 4. Review statistics
cat filtering_report.txt

# 5. Examine output
head -n 5 pathogenic_variants.maf | column -t -s $'\t'
```

## Future Enhancements

Potential future improvements:
- [ ] Julia FFI bindings for calling Rust implementation from Julia
- [ ] Additional predictive scores (CADD, MetaSVM, etc.)
- [ ] VCF output format support
- [ ] Streaming JSON parsing for memory-constrained environments
- [ ] Built-in annotation database updates

## Contributing

This is an internal research tool. For questions, bug reports, or feature requests, please contact:

- Project bioinformaticians
- Development team leads

## Citation

If you use this tool in your research, please cite:

```
JSON2MAF Rust Implementation
Version 0.4.0
Internal Oncology Research Tool
```

## License

Internal use only for oncology research and pharmaceutical development.

## Acknowledgments

- Original Julia implementation by the JSON2MAF development team
- Rust community for excellent libraries:
  - `serde` and `serde_json` for JSON parsing
  - `rayon` for parallel processing
  - `clap` for CLI argument parsing
  - `flate2` for gzip compression
  - `csv` for MAF output
  - `anyhow` and `thiserror` for error handling

---

**Version**: 0.4.0
**Last Updated**: February 2026
**Status**: Production-ready
**Test Coverage**: 25 passing tests
**Supported Rust**: 1.70+
