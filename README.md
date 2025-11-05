# JSON2MAF

A high-performance Julia tool for filtering pathogenic variants from Illumina Nirvana annotation output and converting them to MAF (Mutation Annotation Format).

[![Julia](https://img.shields.io/badge/Julia-1.10+-blue.svg)](https://julialang.org/)
[![License](https://img.shields.io/badge/License-Internal-red.svg)]()
[![Tests](https://img.shields.io/badge/Tests-435%2F435%20passing-brightgreen.svg)]()

## Overview

JSON2MAF processes variant annotation data from Illumina Nirvana (gzipped JSON format) and extracts pathogenic/likely pathogenic variants for oncology research. The tool implements intelligent filtering logic combining ClinVar annotations with multiple predictive scoring systems.

### Key Features

- **Multi-threaded parallel processing** - Utilizes multiple CPU cores for optimal performance
- **Quality filtering** - Filters variants based on sequencing depth and variant allele frequency (VAF)
- **Population frequency filtering** - Excludes common variants in East Asian populations
- **ClinVar prioritization** - Intelligent handling of multiple ClinVar entries with conflict resolution
- **Predictive score integration** - Combines REVEL, DANN, PrimateAI-3D, and COSMIC scores
- **Complete statistics tracking** - Provides detailed reports on filtering stages
- **Configurable thresholds** - All filtering parameters can be customized

### Filtering Decision Logic

1. **ClinVar-based filtering** (Primary):
   - Pathogenic → Include
   - Likely pathogenic → Include
   - Prioritizes higher review status and cancer-related annotations

2. **Predictive score-based filtering** (Secondary):
   - 2+ supporting scores (REVEL, DANN, PrimateAI-3D) → Include as Likely pathogenic
   - PrimateAI-3D alone → Include (highest priority predictive score)

## System Requirements

- **Julia**: 1.10 or higher
- **RAM**: 2-4 GB (varies with input file size and thread count)
- **CPU**: Multi-core processor recommended for optimal performance

### Dependencies

All dependencies are automatically installed via Julia's package manager:

```julia
JSON3.jl          # JSON parsing
CodecZlib.jl      # Gzip file handling
DataFrames.jl     # Data manipulation
CSV.jl            # MAF output
ArgParse.jl       # Command-line interface
ProgressMeter.jl  # Progress display
```

## Installation

```bash
# Clone the repository
git clone https://github.com/your-org/JSON2MAF.git
cd JSON2MAF

# Install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run tests (optional)
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Usage

### Quick Start

```bash
# Basic usage
julia --project=. -t 6 bin/json2maf.jl \
  -i input.json.gz \
  -o output.maf

# With verbose output and statistics report
julia --project=. -t 6 bin/json2maf.jl \
  -i input.json.gz \
  -o output.maf \
  --stats report.txt \
  --verbose
```

### Complete Example with Custom Parameters

```bash
julia --project=. -t 6 bin/json2maf.jl \
  --input input.json.gz \
  --output output.maf \
  --min-depth 30 \
  --min-vaf 0.03 \
  --max-eas-af 0.01 \
  --min-revel 0.75 \
  --min-primate-ai 0.8 \
  --min-dann 0.96 \
  --stats report.txt \
  --verbose
```

### Command-Line Options

**Required Arguments**:
- `-i, --input FILE`: Input Nirvana JSON.gz file path
- `-o, --output FILE`: Output MAF file path

**Quality Filtering Parameters**:
- `--min-depth INT`: Minimum sequencing depth (default: 30)
- `--min-vaf FLOAT`: Minimum variant allele frequency (default: 0.03)
- `--max-eas-af FLOAT`: Maximum East Asian allele frequency (default: 0.01)

**Predictive Score Thresholds**:
- `--min-revel FLOAT`: REVEL score threshold (default: 0.75)
- `--min-primate-ai FLOAT`: PrimateAI-3D score threshold (default: 0.8)
- `--min-dann FLOAT`: DANN score threshold (default: 0.96)

**Other Options**:
- `-v, --verbose`: Enable verbose output
- `-q, --quiet`: Suppress progress bar
- `--stats FILE`: Save statistics report to file
- `--keep-temp`: Keep temporary thread files (for debugging)
- `-h, --help`: Show help message
- `--version`: Show version information

### Performance Tuning

Adjust thread count based on your system:

```bash
# Set thread count before running Julia
export JULIA_NUM_THREADS=8

# Or specify inline
julia --project=. -t 8 bin/json2maf.jl -i input.json.gz -o output.maf
```

**Recommended thread counts**:
- 4 threads: Standard workstations
- 6-8 threads: High-performance workstations
- 16+ threads: Server environments

## Output Format

### MAF File

The tool generates a standard MAF (Mutation Annotation Format) file with the following key fields:

| Field | Description |
|-------|-------------|
| `Hugo_Symbol` | Gene symbol |
| `Chromosome` | Chromosome name |
| `Start_Position` | Variant start position |
| `End_Position` | Variant end position |
| `Variant_Classification` | Consequence type (e.g., Missense_Mutation) |
| `Variant_Type` | Variant type (SNP, DNP, INS, DEL) |
| `Reference_Allele` | Reference allele |
| `Tumor_Seq_Allele2` | Alternate allele |
| `dbSNP_RS` | dbSNP ID |
| `COSMIC_ID` | COSMIC ID |
| `HGVSc` | HGVS coding notation |
| `HGVSp_Short` | HGVS protein notation |
| `Transcript_ID` | Transcript identifier |
| `ClinVar_Significance` | ClinVar clinical significance |
| `ClinVar_Review_Status` | ClinVar review status |
| `ClinVar_Disease` | Associated diseases |
| `REVEL_Score` | REVEL pathogenicity score |
| `DANN_Score` | DANN score |
| `PrimateAI_Score` | PrimateAI-3D score |
| `COSMIC_Count` | Number of COSMIC samples |

### Statistics Report

When using `--stats` or `--verbose`, the tool generates a detailed report:

```
═══════════════════════════════════════════════════════════
                    Filtering Statistics
═══════════════════════════════════════════════════════════

Processing Mode:        Multi-threaded parallel processing
Thread Count:           6

Quality Filtering:
  - Passed quality:     1,340
  - Insufficient depth: 623
  - Low VAF:            1,193
  - High population AF: 27,545

Pathogenicity Assessment:
  - ClinVar Pathogenic:         0
  - ClinVar Likely pathogenic:  0
  - Predictive score support:   198
    * PrimateAI-3D only:        7
    * 2+ scores support:        191

Final Results:
  - Included variants:  198
  - Excluded variants:  1,142

═══════════════════════════════════════════════════════════
```

## Input Requirements

### Nirvana JSON Format

The input must be a gzipped JSON file produced by Illumina Nirvana annotation pipeline. Required fields include:

- `header`: Metadata and annotation sources
- `positions`: Per-variant annotations (main processing target)
  - Basic variant information (chromosome, position, alleles)
  - `samples`: Sample-level data (depth, VAF)
  - `transcripts`: Transcript consequences
  - `clinvar`: ClinVar annotations
  - `populationFrequencies`: gnomAD, 1000 Genomes data
  - Predictive scores: `primateAI-3D`, `DANN`, `REVEL`

### Data Quality Considerations

- Variants must pass VCF FILTER field (only "PASS" variants are processed)
- Missing annotations are handled gracefully
- Boolean fields in Nirvana JSON only appear when true

## Limitations and Constraints

### Technical Limitations

1. **Memory Usage**:
   - Proportional to input file size and thread count
   - Typical usage: 2-4 GB for 30K-50K variants with 6 threads
   - For very large files (>100MB), monitor RAM usage

2. **Thread Scalability**:
   - Optimal speedup with 4-8 threads
   - Diminishing returns beyond 16 threads due to overhead
   - Parallel efficiency: ~87% with 4-6 threads

3. **Input File Format**:
   - Only supports Nirvana JSON format
   - Must be gzipped (`.json.gz` extension)
   - Cannot process VCF files directly (requires Nirvana annotation first)

### Biological/Clinical Limitations

1. **Pathogenicity Assessment**:
   - Limited to ClinVar and computational predictions
   - Does not include functional validation
   - Cancer-relevance depends on ClinVar disease annotations

2. **Predictive Score Availability**:
   - Not all variants have predictive scores
   - Missing scores are treated as non-supporting evidence
   - Score thresholds based on literature recommendations

3. **Population-specific Filtering**:
   - Focused on East Asian populations (EAS AF)
   - Other populations may require parameter adjustments
   - Missing population frequency data is conservatively allowed

4. **Variant Types**:
   - Optimized for SNVs and small indels
   - Large structural variants may not be properly assessed
   - Copy number variations not supported

### Recommended Use Cases

✅ **Appropriate for**:
- Somatic variant calling from tumor samples
- Germline variant analysis for cancer predisposition
- Research cohort variant prioritization
- Clinical biomarker discovery

❌ **Not recommended for**:
- Direct clinical diagnostic reporting (requires expert review)
- Pharmacogenomics applications
- Population genetics studies
- Real-time clinical decision support

### Important Notes

- This tool is designed for **research purposes in oncology drug development**
- Results should be reviewed by qualified bioinformaticians and clinicians
- Filtering parameters may need adjustment based on specific research goals
- ClinVar annotations are periodically updated; consider re-annotation for critical studies

## Performance Benchmarks

Tested on a workstation with Intel Core i7 (4 cores) and 16 GB RAM:

| File Size | Variants | Threads | Processing Time | Memory | Parallel Efficiency |
|-----------|----------|---------|-----------------|--------|-------------------|
| 42 MB | 30,701 | 1 | ~35 min | 1.8 GB | - |
| 42 MB | 30,701 | 4 | 9.8 min | 2.4 GB | 87% |
| 42 MB | 30,701 | 6 | 7.5 min | 2.8 GB | 79% |

*Efficiency = (Sequential Time / Parallel Time) / Thread Count × 100%*

## Project Structure

```
JSON2MAF/
├── src/
│   ├── JSON2MAF.jl              # Main module
│   ├── parser/
│   │   └── NirvanaParser.jl     # JSON parsing
│   ├── filters/
│   │   ├── QualityFilter.jl     # Quality filtering
│   │   ├── ClinVarFilter.jl     # ClinVar filtering
│   │   ├── PredictiveScores.jl  # Predictive scores
│   │   └── DecisionEngine.jl    # Decision logic
│   ├── converters/
│   │   └── MAFConverter.jl      # MAF conversion
│   ├── writers/
│   │   └── MAFWriter.jl         # MAF output
│   └── utils/
│       ├── DataStructures.jl    # Data structures
│       └── FilterConfig.jl      # Configuration
├── bin/
│   └── json2maf.jl              # Main executable
├── test/
│   ├── runtests.jl              # Test runner
│   └── *_tests.jl               # Test suites (435 tests)
└── docs/
    ├── CLAUDE.md                # Development guide
    ├── Requirement.md           # Requirements specification
    └── nirvana_json_format_spec.md  # Format documentation
```

## Testing

Run the comprehensive test suite:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Current status: ✅ **435/435 tests passing**

Test coverage includes:
- Data structure validation
- JSON parsing accuracy
- Quality filtering logic
- ClinVar prioritization
- Predictive score assessment
- Decision engine rules
- MAF format compliance

## Troubleshooting

### Common Issues

**Issue**: Out of memory errors
- **Solution**: Reduce thread count or process smaller batches

**Issue**: Slow processing with multiple threads
- **Solution**: Verify `JULIA_NUM_THREADS` is set correctly

**Issue**: No variants in output
- **Solution**: Check filtering parameters (may be too stringent)

**Issue**: Missing predictive scores
- **Solution**: Ensure Nirvana annotation included all databases

## Contributing

This is an internal research tool. For questions or suggestions, please contact:
- Project bioinformaticians
- Product managers

## Citation

If you use this tool in your research, please cite:
```
JSON2MAF: A high-performance variant filtering tool for oncology research
[Your Organization]
Version 0.4.0 (2024)
```

## License

Internal use only for oncology research and drug development.

## Acknowledgments

- Illumina Nirvana annotation pipeline
- ClinVar database (NCBI)
- gnomAD project
- COSMIC database (Wellcome Sanger Institute)

---

**Version**: 0.4.0
**Last Updated**: November 2024
**Status**: Production-ready
