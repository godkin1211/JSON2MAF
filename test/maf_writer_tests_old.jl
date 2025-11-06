"""
maf_writer_tests.jl

Test MAF file writing functionality
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
        # Create a complete MAFRecord
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
        # Create a MAFRecord with empty values
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

    @testset "Write Single Record" begin
        # Create a test record
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

        # Write to a temporary file
        tmpfile = tempname() * ".maf"
        write_maf_file(tmpfile, record)

        @test isfile(tmpfile)

        # Read and validate (force all columns to be strings)
        df = CSV.read(tmpfile, DataFrame, delim='\t', types=String)
        @test nrow(df) == 1
        @test ncol(df) == 29
        @test df[1, :Hugo_Symbol] == "TP53"
        @test df[1, :Chromosome] == "chr17"
        @test df[1, :Start_Position] == "7579472"
        @test df[1, :Variant_Classification] == "Missense_Mutation"
        @test df[1, :ClinVar_ID] == "RCV000012345"
        @test df[1, :PrimateAI_Score] == "0.85"

        # Clean up
        rm(tmpfile)
    end

    @testset "Write Multiple Records" begin
        # Create multiple test records
        records = [
            MAFRecord(
                hugo_symbol = "TP53",
                chromosome = "chr17",
                start_position = 7579472,
                end_position = 7579472,
                variant_classification = "Missense_Mutation",
                variant_type = "SNP",
                reference_allele = "G",
                tumor_seq_allele2 = "C"
            ),
            MAFRecord(
                hugo_symbol = "BRCA1",
                chromosome = "chr17",
                start_position = 43045677,
                end_position = 43045677,
                variant_classification = "Nonsense_Mutation",
                variant_type = "SNP",
                reference_allele = "C",
                tumor_seq_allele2 = "T"
            ),
            MAFRecord(
                hugo_symbol = "EGFR",
                chromosome = "chr7",
                start_position = 55249071,
                end_position = 55249071,
                variant_classification = "Missense_Mutation",
                variant_type = "SNP",
                reference_allele = "T",
                tumor_seq_allele2 = "G"
            )
        ]

        # Write to a temporary file
        tmpfile = tempname() * ".maf"
        write_maf_file(tmpfile, records)

        @test isfile(tmpfile)

        # Read and validate (force all columns to be strings)
        df = CSV.read(tmpfile, DataFrame, delim='\t', types=String)
        @test nrow(df) == 3
        @test df[1, :Hugo_Symbol] == "TP53"
        @test df[2, :Hugo_Symbol] == "BRCA1"
        @test df[3, :Hugo_Symbol] == "EGFR"
        @test df[2, :Variant_Classification] == "Nonsense_Mutation"

        # Clean up
        rm(tmpfile)
    end

    @testset "Write Empty Record List" begin
        tmpfile = tempname() * ".maf"
        write_maf_file(tmpfile, MAFRecord[])

        @test isfile(tmpfile)

        # Should only have the header
        lines = readlines(tmpfile)
        @test length(lines) == 1
        @test startswith(lines[1], "Hugo_Symbol\t")

        # Clean up
        rm(tmpfile)
    end

    @testset "Append Mode" begin
        tmpfile = tempname() * ".maf"

        # Write the first record
        record1 = MAFRecord(
            hugo_symbol = "TP53",
            chromosome = "chr17",
            start_position = 100,
            end_position = 100
        )
        write_maf_file(tmpfile, record1)

        # Append the second record
        record2 = MAFRecord(
            hugo_symbol = "BRCA1",
            chromosome = "chr17",
            start_position = 200,
            end_position = 200
        )
        write_maf_file(tmpfile, record2, append=true)

        # Validate that there are two records (force all columns to be strings)
        df = CSV.read(tmpfile, DataFrame, delim='\t', types=String)
        @test nrow(df) == 2
        @test df[1, :Hugo_Symbol] == "TP53"
        @test df[2, :Hugo_Symbol] == "BRCA1"

        # Clean up
        rm(tmpfile)
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

        tmpfile = tempname() * ".maf"
        write_maf_file(tmpfile, record)

        # Read the first line (header)
        lines = readlines(tmpfile)
        header = split(lines[1], '\t')

        # Validate the position of important columns
        @test header[1] == "Hugo_Symbol"
        @test header[23] == "PrimateAI_Score"
        @test header[24] == "DANN_Score"
        @test header[25] == "REVEL_Score"
        @test header[26] == "gnomAD_AF"
        @test header[29] == "VAF"

        # Clean up
        rm(tmpfile)
    end
end
