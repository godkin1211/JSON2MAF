"""
maf_converter_tests.jl

Test MAF format conversion functionality
"""

using Test
using JSON2MAF

@testset "MAF Format Conversion Tests" begin

    @testset "Variant Classification Mapping" begin
        @test map_variant_classification("missense_variant") == "Missense_Mutation"
        @test map_variant_classification("stop_gained") == "Nonsense_Mutation"
        @test map_variant_classification("frameshift_variant") == "Frame_Shift_Indel"
        @test map_variant_classification("splice_donor_variant") == "Splice_Site"
        @test map_variant_classification("splice_acceptor_variant") == "Splice_Site"
        @test map_variant_classification("synonymous_variant") == "Silent"
        @test map_variant_classification("3_prime_UTR_variant") == "3'UTR"
        @test map_variant_classification("5_prime_UTR_variant") == "5'UTR"
        @test map_variant_classification("intron_variant") == "Intron"
        @test map_variant_classification("intergenic_variant") == "IGR"
        @test map_variant_classification("start_lost") == "Translation_Start_Site"
        @test map_variant_classification("stop_lost") == "Nonstop_Mutation"
        @test map_variant_classification("inframe_deletion") == "In_Frame_Del"
        @test map_variant_classification("inframe_insertion") == "In_Frame_Ins"
    end

    @testset "Variant Type Mapping" begin
        # SNP
        @test map_variant_type("A", "G") == "SNP"
        @test map_variant_type("C", "T") == "SNP"

        # DNP
        @test map_variant_type("AA", "GG") == "DNP"
        @test map_variant_type("CT", "AG") == "DNP"

        # TNP
        @test map_variant_type("AAA", "GGG") == "TNP"

        # ONP
        @test map_variant_type("AAAA", "GGGG") == "ONP"

        # Insertion
        @test map_variant_type("A", "ATG") == "INS"
        @test map_variant_type("", "AAA") == "INS"

        # Deletion
        @test map_variant_type("ATG", "A") == "DEL"
        @test map_variant_type("AAA", "") == "DEL"
    end

    @testset "HGVS Extraction" begin
        transcript = TranscriptAnnotation(
            "ENST00000123456",  # id
            "TP53",              # gene_symbol
            nothing,             # hgnc
            ["missense_variant"], # consequence
            "V/M",               # amino_acids
            100,                 # cdna_pos
            99,                  # cds_pos
            33,                  # protein_pos
            "c.215C>G",          # hgvsc
            "p.Pro72Arg"         # hgvsp
        )

        hgvsc, hgvsp, hgvsp_short = extract_hgvs_notation(transcript)

        @test hgvsc == "c.215C>G"
        @test hgvsp == "p.Pro72Arg"
        @test hgvsp_short != "."  # Should have a value
    end

    @testset "Canonical Transcript Selection" begin
        # Test empty list
        @test select_canonical_transcript(TranscriptAnnotation[]) === nothing

        # Test single transcript
        t1 = TranscriptAnnotation(
            "ENST00000123456", "TP53", nothing, ["missense_variant"],
            "V/M", 100, 99, 33, "c.215C>G", "p.Pro72Arg"
        )
        @test select_canonical_transcript([t1]) === t1

        # Test multiple transcripts - Nirvana places canonical first
        t2 = TranscriptAnnotation(
            "ENST00000222222", "TP53", nothing, ["missense_variant"],
            "V/M", 100, 99, 33, "c.215C>G", "p.Pro72Arg"
        )
        t3 = TranscriptAnnotation(
            "ENST00000333333", "TP53", nothing, ["missense_variant"],
            "V/M", 100, 99, 33, "c.215C>G", "p.Pro72Arg"
        )

        # Should select the first one (canonical)
        @test select_canonical_transcript([t2, t1, t3]) === t2

        # Always select the first one
        @test select_canonical_transcript([t1, t3]) === t1
    end

    @testset "Complete MAF Conversion" begin
        # Create test variant
        transcript = TranscriptAnnotation(
            "ENST00000269305",  # id
            "TP53",              # gene_symbol
            nothing,             # hgnc
            ["missense_variant"], # consequence
            "V/M",               # amino_acids
            215,                 # cdna_pos
            215,                 # cds_pos
            72,                  # protein_pos
            "c.215C>G",          # hgvsc
            "p.Pro72Arg"         # hgvsp
        )

        clinvar_entry = ClinVarEntry(
            "RCV000012345",
            "12345",
            "Pathogenic",
            "reviewed by expert panel",
            ["Li-Fraumeni syndrome"],
            "2024-01-01"
        )

        cosmic_entry = CosmicEntry(
            "COSM43598",
            "TP53",
            "missense",
            150
        )

        pop_freq = PopulationFrequency(
            "gnomad-exome",
            0.02,   # all_af
            0.01,   # eas_af
            nothing, nothing, nothing
        )

        variant = VariantPosition(
            "chr17",
            7579472,
            7579472,
            "G",
            "C",
            "SNV",
            150,      # total_depth
            [0.45],   # variant_frequencies
            [transcript],
            [clinvar_entry],
            [cosmic_entry],
            [pop_freq],
            0.85,     # primate_ai_3d
            nothing,
            0.98,     # dann_score
            0.80,     # revel_score
            ["rs1042522"]
        )

        decision = FilterDecision(
            true,
            "Pathogenic",
            "ClinVar",
            "ClinVar: Pathogenic"
        )

        maf_record = variant_to_maf(variant, decision)

        # Verify basic fields
        @test maf_record.hugo_symbol == "TP53"
        @test maf_record.chromosome == "chr17"
        @test maf_record.start_position == 7579472
        @test maf_record.end_position == 7579472
        @test maf_record.reference_allele == "G"
        @test maf_record.tumor_seq_allele2 == "C"

        # Verify variant classification
        @test maf_record.variant_classification == "Missense_Mutation"
        @test maf_record.variant_type == "SNP"

        # Verify HGVS
        @test maf_record.hgvsc == "c.215C>G"
        @test maf_record.hgvsp == "p.Pro72Arg"

        # Verify Transcript ID
        @test maf_record.transcript_id == "ENST00000269305"

        # Verify ClinVar
        @test maf_record.clinvar_id == "RCV000012345"
        @test maf_record.clinvar_review_status == "reviewed by expert panel"
        @test maf_record.clinvar_significance == "Pathogenic"
        @test occursin("Li-Fraumeni", maf_record.clinvar_disease)

        # Verify COSMIC
        @test maf_record.cosmic_id == "COSM43598"

        # Verify dbSNP
        @test maf_record.dbsnp_rs == "rs1042522"

        # Verify predictive scores (should be strings)
        @test maf_record.primate_ai_score == "0.85"
        @test maf_record.dann_score == "0.98"
        @test maf_record.revel_score == "0.8"  # Julia string() removes trailing 0

        # Verify gnomAD AF
        @test maf_record.gnomad_af == "0.02"
        @test maf_record.gnomad_eas_af == "0.01"

        # Verify depth and VAF
        @test maf_record.depth == "150"
        @test maf_record.vaf == "0.45"
    end

    @testset "Missing Field Handling" begin
        # Create variant without transcripts
        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[],  # Empty transcripts
            ClinVarEntry[],
            CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )

        decision = FilterDecision(true, "Likely pathogenic", "Predictive", "Reason")

        maf_record = variant_to_maf(variant, decision)

        # Verify missing values are handled as "."
        @test maf_record.hugo_symbol == "."
        @test maf_record.variant_classification == "."
        @test maf_record.hgvsc == "."
        @test maf_record.hgvsp == "."
        @test maf_record.transcript_id == "."
        @test maf_record.clinvar_id == "."
        @test maf_record.cosmic_id == "."
        @test maf_record.dbsnp_rs == "."
        @test maf_record.primate_ai_score === nothing
        @test maf_record.gnomad_af === nothing
    end

    @testset "Multiple COSMIC ID Handling" begin
        cosmic1 = CosmicEntry("COSM111", "TP53", "missense", 50)
        cosmic2 = CosmicEntry("COSM222", "TP53", "missense", 30)

        variant = VariantPosition(
            "chr17", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[],
            ClinVarEntry[],
            [cosmic1, cosmic2],  # Multiple COSMIC
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )

        decision = FilterDecision(true, "Likely pathogenic", "Predictive", "Reason")
        maf_record = variant_to_maf(variant, decision)

        # Should separate multiple IDs with ;
        @test occursin("COSM111", maf_record.cosmic_id)
        @test occursin("COSM222", maf_record.cosmic_id)
        @test occursin(";", maf_record.cosmic_id)
    end

    @testset "Multiple Disease Name Handling" begin
        clinvar_entry = ClinVarEntry(
            "RCV12345", "12345",
            "Pathogenic",
            "criteria provided, multiple submitters, no conflicts",
            ["Disease A", "Disease B", "Disease C"],
            "2024-01-01"
        )

        variant = VariantPosition(
            "chr1", 100, 100, "A", "G", "SNV",
            100, [0.5],
            TranscriptAnnotation[],
            [clinvar_entry],
            CosmicEntry[],
            PopulationFrequency[],
            nothing, nothing, nothing, nothing,
            String[]
        )

        decision = FilterDecision(true, "Pathogenic", "ClinVar", "Reason")
        maf_record = variant_to_maf(variant, decision)

        # Disease names should be separated with ;
        @test occursin("Disease A", maf_record.clinvar_disease)
        @test occursin("Disease B", maf_record.clinvar_disease)
        @test occursin(";", maf_record.clinvar_disease)
    end
end
