"""
maf_writer_tests.jl

Test MAF file writing functionality

NOTE: The write_maf_file() function and format_contributing_scores() function have been removed.
The current API uses batch writing functions: create_maf_writer, write_maf_batch,
close_maf_writer, and merge_maf_files.
"""

using Test
using JSON2MAF
using DataFrames
using CSV

@testset "MAF Writer Tests" begin

    @testset "MAF Header" begin
        header = maf_header()

        @test length(header) == 29
        @test header[1] == "Hugo_Symbol"
        @test header[2] == "Chromosome"
        @test header[end] == "VAF"
        @test "ClinVar_ID" in header
        @test "PrimateAI_Score" in header
    end

    @testset "MAFRecord to Row Conversion" begin
        # Create complete MAFRecord
        record = MAFRecord(
            hugo_symbol = "TP53",
            chromosome = "chr17",
            start_position = 7579472,
            end_position = 7579472,
            strand = "+",
            variant_classification = "Missense_Mutation",
            variant_type = "SNP",
            reference_allele = "G",
            tumor_seq_allele1 = "G",
            tumor_seq_allele2 = "C",
            tumor_sample_barcode = "SAMPLE001",
            hgvsc = "c.215C>G",
            hgvsp = "p.Pro72Arg",
            hgvsp_short = "p.P72R",
            transcript_id = "ENST00000269305",
            dbsnp_rs = "rs1042522",
            dbsnp_val_status = "",
            cosmic_id = "COSM43598",
            clinvar_id = "RCV000012345",
            clinvar_review_status = "reviewed by expert panel",
            clinvar_significance = "Pathogenic",
            clinvar_disease = "Li-Fraumeni syndrome",
            primate_ai_score = "0.85",
            dann_score = "0.98",
            revel_score = "0.8",
            gnomad_af = "0.02",
            gnomad_eas_af = "0.01",
            depth = "150",
            vaf = "0.45"
        )

        row = mafrecord_to_row(record)

        @test length(row) == 29
        @test row[1] == "TP53"
        @test row[2] == "chr17"
        @test row[3] == "7579472"
        @test row[6] == "Missense_Mutation"
        @test row[12] == "c.215C>G"
        @test row[19] == "RCV000012345"
        @test row[23] == "0.85"
    end

    @testset "Empty Value Formatting" begin
        # Create MAFRecord with empty values
        record = MAFRecord(
            hugo_symbol = "TP53",
            chromosome = "chr17",
            start_position = 100,
            end_position = 100,
            strand = "+",
            variant_classification = "Missense_Mutation",
            variant_type = "SNP",
            reference_allele = "A",
            tumor_seq_allele1 = "A",
            tumor_seq_allele2 = "G",
            tumor_sample_barcode = "",
            hgvsc = "c.100A>G",
            hgvsp = "",
            hgvsp_short = "",
            transcript_id = "ENST123",
            dbsnp_rs = "",
            dbsnp_val_status = "",
            cosmic_id = "",
            clinvar_id = "",
            clinvar_review_status = "",
            clinvar_significance = "",
            clinvar_disease = "",
            primate_ai_score = nothing,
            dann_score = nothing,
            revel_score = nothing,
            gnomad_af = nothing,
            gnomad_eas_af = nothing,
            depth = "",
            vaf = ""
        )

        row = mafrecord_to_row(record)

        # Empty values should be converted to "."
        @test row[11] == "."  # tumor_sample_barcode
        @test row[13] == "."  # hgvsp
        @test row[16] == "."  # dbsnp_rs
        @test row[18] == "."  # cosmic_id
        @test row[19] == "."  # clinvar_id
        @test row[23] == "."  # primate_ai_score (nothing)
        @test row[24] == "."  # dann_score (nothing)
        @test row[26] == "."  # gnomad_af (nothing)
    end

    # NOTE: Tests for write_maf_file() function have been removed because this function
    # no longer exists in the codebase. The new API uses batch writing functions:
    # - create_maf_writer() - Creates a writer for batch operations
    # - write_maf_batch() - Writes a batch of records
    # - close_maf_writer() - Closes the writer
    # - merge_maf_files() - Merges multiple MAF files

    @testset "Legacy write_maf_file Tests (Disabled)" begin
        @test_skip "Write single record - write_maf_file() no longer exists, use batch writer API"
        @test_skip "Write multiple records - write_maf_file() no longer exists, use batch writer API"
        @test_skip "Write empty record list - write_maf_file() no longer exists, use batch writer API"
        @test_skip "Append mode - write_maf_file() no longer exists, use batch writer API"
    end

    @testset "Column Order Correctness" begin
        record = MAFRecord(
            hugo_symbol = "TEST",
            chromosome = "chr1",
            start_position = 1,
            end_position = 1,
            primate_ai_score = "0.9",
            dann_score = "0.95",
            revel_score = "0.8"
        )

        row = mafrecord_to_row(record)

        # Verify important field positions by checking header
        header = maf_header()
        @test header[1] == "Hugo_Symbol"
        @test header[23] == "PrimateAI_Score"
        @test header[24] == "DANN_Score"
        @test header[25] == "REVEL_Score"
        @test header[26] == "gnomAD_AF"
        @test header[29] == "VAF"
    end

    # Batch Writer API Tests
    @testset "Batch Writer API" begin
        # Test that the batch writer functions exist
        @test isdefined(JSON2MAF, :create_maf_writer)
        @test isdefined(JSON2MAF, :write_maf_batch)
        @test isdefined(JSON2MAF, :close_maf_writer)
        @test isdefined(JSON2MAF, :merge_maf_files)

        # Create a test record
        record = MAFRecord(
            hugo_symbol = "BRCA1",
            chromosome = "chr17",
            start_position = 43045677,
            end_position = 43045677,
            variant_classification = "Missense_Mutation",
            variant_type = "SNP",
            reference_allele = "G",
            tumor_seq_allele2 = "A"
        )

        # Test basic batch writer workflow
        tmpfile = tempname() * ".maf"

        try
            writer = create_maf_writer(tmpfile)
            @test writer !== nothing

            write_maf_batch(writer, [record])
            close_maf_writer(writer)

            # Verify file was created and has content
            @test isfile(tmpfile)
            df = CSV.read(tmpfile, DataFrame, delim='\t', types=String)
            @test nrow(df) == 1
            @test df[1, :Hugo_Symbol] == "BRCA1"

            # Clean up
            rm(tmpfile)
        catch e
            # Clean up on error
            isfile(tmpfile) && rm(tmpfile)
            rethrow(e)
        end
    end
end
