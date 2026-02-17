/// Integration tests for JSON2MAF
/// Tests end-to-end parsing, filtering, and MAF conversion

use json2maf::*;
use std::fs::File;
use std::io::Write;
use tempfile::TempDir;
use flate2::write::GzEncoder;
use flate2::Compression;

#[test]
fn test_parse_and_convert_with_new_fields() {
    // Test that new annotation fields are correctly parsed and converted to MAF
    let test_json = r#"{
        "header": {
            "annotator": "Nirvana 3.0",
            "creationTime": "2024-01-01",
            "genomeAssembly": "GRCh38",
            "schemaVersion": 6,
            "dataSources": [],
            "samples": ["TEST"]
        },
        "positions": [{
            "chromosome": "chr7",
            "position": 140453136,
            "refAllele": "A",
            "altAlleles": ["T"],
            "filters": ["PASS"],
            "samples": [{
                "totalDepth": 100,
                "variantFrequencies": [0.45]
            }],
            "variants": [{
                "variantType": "SNV",
                "transcripts": [{
                    "transcript": "NM_004333.4",
                    "hgnc": "BRAF",
                    "consequence": ["missense_variant"],
                    "impact": "moderate",
                    "aminoAcids": "V/E",
                    "cdnaPos": "1799/2301",
                    "cdsPos": "1799/2301",
                    "exons": "15/18",
                    "codons": "Gtg/Gag",
                    "proteinPos": "600/766",
                    "hgvsc": "NM_004333.4:c.1799T>A",
                    "hgvsp": "NP_004324.2:p.Val600Glu",
                    "isCanonical": true
                }],
                "clinvar": [{
                    "id": "RCV000012345",
                    "significance": ["Pathogenic"],
                    "reviewStatus": "criteria provided, multiple submitters",
                    "phenotypes": ["Melanoma", "Colorectal cancer"]
                }]
            }]
        }]
    }"#;

    // Create temporary files
    let temp_dir = TempDir::new().unwrap();
    let input_path = temp_dir.path().join("test.json.gz");
    let _output_path = temp_dir.path().join("output.maf");

    // Write gzipped JSON
    let file = File::create(&input_path).unwrap();
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder.write_all(test_json.as_bytes()).unwrap();
    encoder.finish().unwrap();

    // Parse the JSON
    let (header, variants) = parser::parse_nirvana_json(input_path.to_str().unwrap()).unwrap();

    assert_eq!(header.genome_assembly, "GRCh38");
    assert_eq!(variants.len(), 1);

    let variant = &variants[0];

    // Verify transcript fields are parsed correctly
    assert_eq!(variant.transcripts.len(), 1);
    let transcript = &variant.transcripts[0];

    assert_eq!(transcript.hgnc.as_deref(), Some("BRAF"));
    assert_eq!(transcript.impact.as_deref(), Some("moderate"));
    assert_eq!(transcript.amino_acids.as_deref(), Some("V/E"));
    assert_eq!(transcript.cdna_pos.as_deref(), Some("1799/2301"));
    assert_eq!(transcript.cds_pos.as_deref(), Some("1799/2301"));
    assert_eq!(transcript.exons.as_deref(), Some("15/18"));
    assert_eq!(transcript.codons.as_deref(), Some("Gtg/Gag"));
    assert_eq!(transcript.protein_pos.as_deref(), Some("600/766"));

    // Test filtering and conversion
    let config = FilterConfig::default();

    let quality_result = filters::quality::apply_quality_filters(variant, &config);
    assert!(quality_result.passes_quality);

    let clinvar_result = filters::clinvar::assess_clinvar_pathogenicity(&variant.clinvar);
    assert!(clinvar_result.is_pathogenic);

    let predictive_result = filters::predictive::assess_predictive_scores(variant, &config);

    let decision = filters::decision::make_filter_decision(
        variant,
        &clinvar_result,
        &predictive_result,
    );
    assert!(decision.should_include);
    assert_eq!(decision.pathogenicity_class, "Pathogenic");

    // Convert to MAF and verify new fields
    let maf_record = converter::variant_to_maf(variant, &decision);

    assert_eq!(maf_record.hugo_symbol, "BRAF");
    assert_eq!(maf_record.chromosome, "chr7");
    assert_eq!(maf_record.variant_classification, "Missense_Mutation");

    // Verify new fields
    assert_eq!(maf_record.impact, "MODERATE"); // Should be uppercase
    assert_eq!(maf_record.consequence, "missense_variant");
    assert_eq!(maf_record.exon, "15/18");
    assert_eq!(maf_record.codons, "Gtg/Gag");
    assert_eq!(maf_record.amino_acids, "V/E");
    assert_eq!(maf_record.cdna_position, "1799/2301");
    assert_eq!(maf_record.cds_position, "1799/2301");
    assert_eq!(maf_record.protein_position, "600/766");

    // Verify HGVS fields
    assert_eq!(maf_record.hgvsc, "NM_004333.4:c.1799T>A");
    assert_eq!(maf_record.hgvsp, "NP_004324.2:p.Val600Glu");
    assert_eq!(maf_record.hgvsp_short, "p.V600E");

    // Verify ClinVar fields
    assert_eq!(maf_record.clinvar_id, "RCV000012345");
    assert_eq!(maf_record.clinvar_significance, "Pathogenic");
}

#[test]
fn test_missing_annotation_fields() {
    // Test that missing fields are handled gracefully
    let test_json = r#"{
        "header": {
            "annotator": "Nirvana 3.0",
            "creationTime": "2024-01-01",
            "genomeAssembly": "GRCh38",
            "schemaVersion": 6,
            "dataSources": [],
            "samples": ["TEST"]
        },
        "positions": [{
            "chromosome": "chr1",
            "position": 12345,
            "refAllele": "A",
            "altAlleles": ["G"],
            "filters": ["PASS"],
            "samples": [{
                "totalDepth": 50,
                "variantFrequencies": [0.3]
            }],
            "variants": [{
                "variantType": "SNV",
                "transcripts": [{
                    "transcript": "NM_001234.1",
                    "hgnc": "GENE1",
                    "consequence": ["synonymous_variant"]
                }]
            }]
        }]
    }"#;

    let temp_dir = TempDir::new().unwrap();
    let input_path = temp_dir.path().join("test.json.gz");
    let file = File::create(&input_path).unwrap();
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder.write_all(test_json.as_bytes()).unwrap();
    encoder.finish().unwrap();

    let (_, variants) = parser::parse_nirvana_json(input_path.to_str().unwrap()).unwrap();
    assert_eq!(variants.len(), 1);

    let variant = &variants[0];
    let transcript = &variant.transcripts[0];

    // Verify missing fields are None
    assert!(transcript.impact.is_none());
    assert!(transcript.amino_acids.is_none());
    assert!(transcript.cdna_pos.is_none());
    assert!(transcript.exons.is_none());
    assert!(transcript.codons.is_none());

    // Convert to MAF - should use empty strings for missing fields
    let config = FilterConfig::default();
    let clinvar = filters::clinvar::assess_clinvar_pathogenicity(&variant.clinvar);
    let predictive = filters::predictive::assess_predictive_scores(variant, &config);
    let decision = filters::decision::make_filter_decision(variant, &clinvar, &predictive);

    let maf_record = converter::variant_to_maf(variant, &decision);

    assert_eq!(maf_record.impact, "");
    assert_eq!(maf_record.amino_acids, "");
    assert_eq!(maf_record.cdna_position, "");
    assert_eq!(maf_record.exon, "");
    assert_eq!(maf_record.codons, "");
}

#[test]
fn test_multiple_consequences() {
    // Test that multiple consequences are joined with commas
    let test_json = r#"{
        "header": {
            "annotator": "Nirvana 3.0",
            "creationTime": "2024-01-01",
            "genomeAssembly": "GRCh38",
            "schemaVersion": 6,
            "dataSources": [],
            "samples": ["TEST"]
        },
        "positions": [{
            "chromosome": "chr1",
            "position": 12345,
            "refAllele": "A",
            "altAlleles": ["G"],
            "filters": ["PASS"],
            "samples": [{
                "totalDepth": 50,
                "variantFrequencies": [0.3]
            }],
            "variants": [{
                "variantType": "SNV",
                "transcripts": [{
                    "transcript": "NM_001234.1",
                    "hgnc": "GENE1",
                    "consequence": ["missense_variant", "splice_region_variant"],
                    "impact": "moderate"
                }]
            }]
        }]
    }"#;

    let temp_dir = TempDir::new().unwrap();
    let input_path = temp_dir.path().join("test.json.gz");
    let file = File::create(&input_path).unwrap();
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder.write_all(test_json.as_bytes()).unwrap();
    encoder.finish().unwrap();

    let (_, variants) = parser::parse_nirvana_json(input_path.to_str().unwrap()).unwrap();
    let variant = &variants[0];

    let config = FilterConfig::default();
    let clinvar = filters::clinvar::assess_clinvar_pathogenicity(&variant.clinvar);
    let predictive = filters::predictive::assess_predictive_scores(variant, &config);
    let decision = filters::decision::make_filter_decision(variant, &clinvar, &predictive);

    let maf_record = converter::variant_to_maf(variant, &decision);

    // Multiple consequences should be joined with comma
    assert_eq!(maf_record.consequence, "missense_variant,splice_region_variant");

    // Variant classification should use the first (most severe) consequence
    assert_eq!(maf_record.variant_classification, "Missense_Mutation");
}

#[test]
fn test_impact_case_conversion() {
    // Test that impact is converted to uppercase
    let impacts = vec![
        ("low", "LOW"),
        ("moderate", "MODERATE"),
        ("high", "HIGH"),
        ("modifier", "MODIFIER"),
    ];

    for (i, (input_impact, expected_output)) in impacts.iter().enumerate() {
        let test_json = format!(r#"{{
            "header": {{
                "annotator": "Nirvana 3.0",
                "creationTime": "2024-01-01",
                "genomeAssembly": "GRCh38",
                "schemaVersion": 6,
                "dataSources": [],
                "samples": ["TEST"]
            }},
            "positions": [{{
                "chromosome": "chr1",
                "position": {},
                "refAllele": "A",
                "altAlleles": ["G"],
                "filters": ["PASS"],
                "samples": [{{"totalDepth": 50, "variantFrequencies": [0.3]}}],
                "variants": [{{
                    "variantType": "SNV",
                    "transcripts": [{{
                        "transcript": "NM_001234.1",
                        "hgnc": "GENE1",
                        "consequence": ["missense_variant"],
                        "impact": "{}"
                    }}]
                }}]
            }}]
        }}"#, 12345 + i, input_impact);  // Unique position for each test

        let temp_dir = TempDir::new().unwrap();
        let input_path = temp_dir.path().join(format!("test_{}.json.gz", i));
        let file = File::create(&input_path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(test_json.as_bytes()).unwrap();
        encoder.finish().unwrap();

        // Add a small delay to ensure file system consistency
        std::thread::sleep(std::time::Duration::from_millis(10));

        let (_, variants) = parser::parse_nirvana_json(input_path.to_str().unwrap()).unwrap();
        let variant = &variants[0];

        let config = FilterConfig::default();
        let clinvar = filters::clinvar::assess_clinvar_pathogenicity(&variant.clinvar);
        let predictive = filters::predictive::assess_predictive_scores(variant, &config);
        let decision = filters::decision::make_filter_decision(variant, &clinvar, &predictive);

        let maf_record = converter::variant_to_maf(variant, &decision);

        assert_eq!(maf_record.impact, *expected_output,
            "Impact '{}' should be converted to '{}'", input_impact, expected_output);
    }
}
