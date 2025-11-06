"""
MAFWriter.jl

MAF file writing module
Writes MAFRecord to standard MAF format file (TSV)
"""

module MAFWriter

using CSV
using DataFrames
using ..DataStructures

export maf_header, mafrecord_to_row, create_maf_writer, write_maf_batch, close_maf_writer, get_total_written, merge_maf_files, count_maf_records

"""
Standard MAF file column order
"""
const MAF_COLUMNS = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Strand",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
    "HGVSc",
    "HGVSp",
    "HGVSp_Short",
    "Transcript_ID",
    "dbSNP_RS",
    "dbSNP_Val_Status",
    "COSMIC_ID",
    "ClinVar_ID",
    "ClinVar_Review_Status",
    "ClinVar_Significance",
    "ClinVar_Disease",
    "PrimateAI_Score",
    "DANN_Score",
    "REVEL_Score",
    "gnomAD_AF",
    "gnomAD_EAS_AF",
    "Depth",
    "VAF"
]

"""
    maf_header() -> Vector{String}

Return standard MAF file header
"""
function maf_header()
    return MAF_COLUMNS
end

"""
    format_maf_value(val) -> String

Format MAF field value
- nothing or empty string → "."
- Other values → Convert to string
"""
function format_maf_value(val)
    if val === nothing || val == ""
        return "."
    else
        return string(val)
    end
end

"""
    mafrecord_to_row(record::MAFRecord) -> Vector{String}

Convert MAFRecord to one row of MAF file data

# Parameters
- `record`: MAFRecord object

# Returns
- `Vector{String}`: Field values arranged in MAF_COLUMNS order
"""
function mafrecord_to_row(record::MAFRecord)::Vector{String}
    return [
        format_maf_value(record.hugo_symbol),
        format_maf_value(record.chromosome),
        format_maf_value(record.start_position),
        format_maf_value(record.end_position),
        format_maf_value(record.strand),
        format_maf_value(record.variant_classification),
        format_maf_value(record.variant_type),
        format_maf_value(record.reference_allele),
        format_maf_value(record.tumor_seq_allele1),
        format_maf_value(record.tumor_seq_allele2),
        format_maf_value(record.tumor_sample_barcode),
        format_maf_value(record.hgvsc),
        format_maf_value(record.hgvsp),
        format_maf_value(record.hgvsp_short),
        format_maf_value(record.transcript_id),
        format_maf_value(record.dbsnp_rs),
        format_maf_value(record.dbsnp_val_status),
        format_maf_value(record.cosmic_id),
        format_maf_value(record.clinvar_id),
        format_maf_value(record.clinvar_review_status),
        format_maf_value(record.clinvar_significance),
        format_maf_value(record.clinvar_disease),
        format_maf_value(record.primate_ai_score),
        format_maf_value(record.dann_score),
        format_maf_value(record.revel_score),
        format_maf_value(record.gnomad_af),
        format_maf_value(record.gnomad_eas_af),
        format_maf_value(record.depth),
        format_maf_value(record.vaf)
    ]
end

# ============================================================================
# Batch write API - for optimizing large file processing performance
# ============================================================================

"""
    MAFBatchWriter

Batch MAF writer for optimizing write performance of large numbers of records

# Fields
- `output_path`: Output file path
- `io`: File handle
- `buffer`: Record buffer
- `batch_size`: Batch size
- `total_written`: Total number of records written
"""
mutable struct MAFBatchWriter
    output_path::String
    io::Union{IO, Nothing}
    buffer::Vector{MAFRecord}
    batch_size::Int
    total_written::Int
end

"""
    create_maf_writer(output_path::String; batch_size::Int=1000) -> MAFBatchWriter

Create batch MAF writer

# Parameters
- `output_path`: Output file path
- `batch_size`: Batch size (default 1000), write after accumulating this many records

# Returns
- `MAFBatchWriter`: Batch writer instance

# Usage
```julia
writer = create_maf_writer("output.maf", batch_size=1000)
for record in records
    write_maf_batch(writer, record)
end
close_maf_writer(writer)  # Don't forget to close!
```

# Performance advantages
- Reduces I/O operations
- Batch processing improves efficiency
- Suitable for handling large numbers of records
"""
function create_maf_writer(output_path::String; batch_size::Int=1000)::MAFBatchWriter
    # Open file and write header
    io = open(output_path, "w")
    println(io, join(maf_header(), "\t"))

    return MAFBatchWriter(
        output_path,
        io,
        MAFRecord[],
        batch_size,
        0
    )
end

"""
    write_maf_batch(writer::MAFBatchWriter, record::MAFRecord)

Add record to batch writer

Automatically writes to file when buffer reaches batch_size.

# Parameters
- `writer`: Batch writer
- `record`: MAF record to write
"""
function write_maf_batch(writer::MAFBatchWriter, record::MAFRecord)
    push!(writer.buffer, record)

    # When buffer is full, write to file
    if length(writer.buffer) >= writer.batch_size
        flush_buffer(writer)
    end
end

"""
    flush_buffer(writer::MAFBatchWriter)

Write buffered records to file
"""
function flush_buffer(writer::MAFBatchWriter)
    if isempty(writer.buffer)
        return
    end

    # Convert all buffered records to rows and write
    for record in writer.buffer
        row = mafrecord_to_row(record)
        println(writer.io, join(row, "\t"))
    end

    writer.total_written += length(writer.buffer)
    empty!(writer.buffer)  # Clear buffer
end

"""
    close_maf_writer(writer::MAFBatchWriter)

Close batch writer and ensure all data is written

# Parameters
- `writer`: Batch writer

# Notes
Must call this function to ensure all buffered records are written to file!
"""
function close_maf_writer(writer::MAFBatchWriter)
    # Write remaining buffered records
    flush_buffer(writer)

    # Close file
    if writer.io !== nothing
        close(writer.io)
        writer.io = nothing
    end
end

"""
    get_total_written(writer::MAFBatchWriter) -> Int

Return total number of records written
"""
function get_total_written(writer::MAFBatchWriter)::Int
    return writer.total_written + length(writer.buffer)
end

# ============================================================================
# MAF file merging - for multi-threaded parallel processing
# ============================================================================

"""
    merge_maf_files(input_files::Vector{String}, output_file::String; keep_temp::Bool=false)

Merge multiple MAF files into a single file

This function is used to merge multiple MAF files generated by multi-threaded parallel processing.
Preserves all records and ensures only one standard header.

# Parameters
- `input_files`: List of input MAF file paths
- `output_file`: Output merged MAF file path
- `keep_temp`: Whether to keep temporary files (default false, will delete)

# Usage example
```julia
# Merge after parallel processing
merge_maf_files(
    ["output_1.maf", "output_2.maf", "output_3.maf"],
    "final_output.maf",
    keep_temp=false
)
```

# Notes
- Input files must all exist
- All input files must have the same field structure
- Default deletes temporary files to save space
"""
function merge_maf_files(input_files::Vector{String}, output_file::String; keep_temp::Bool=false)
    # Check input files
    for file in input_files
        if !isfile(file)
            error("Input file not found: $file")
        end
    end

    # Open output file and write header
    open(output_file, "w") do out_io
        # Write standard header
        println(out_io, join(maf_header(), "\t"))

        # Read and merge each input file one by one
        for (idx, input_file) in enumerate(input_files)
            open(input_file, "r") do in_io
                # Skip input file header
                readline(in_io)

                # Copy all data rows
                while !eof(in_io)
                    line = readline(in_io)
                    if !isempty(strip(line))
                        println(out_io, line)
                    end
                end
            end
        end
    end

    # Delete temporary files (if needed)
    if !keep_temp
        for file in input_files
            try
                rm(file)
            catch e
                @warn "Failed to delete temporary file: $file" exception=e
            end
        end
    end

    return output_file
end

"""
    count_maf_records(filepath::String) -> Int

Count number of records in MAF file (excluding header)

# Parameters
- `filepath`: MAF file path

# Returns
- `Int`: Number of records
"""
function count_maf_records(filepath::String)::Int
    if !isfile(filepath)
        error("File not found: $filepath")
    end

    count = 0
    open(filepath, "r") do io
        # Skip header
        readline(io)

        # Count non-empty lines
        while !eof(io)
            line = readline(io)
            if !isempty(strip(line))
                count += 1
            end
        end
    end

    return count
end

end # module MAFWriter
