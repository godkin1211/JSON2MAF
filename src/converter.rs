use crate::types::*;

pub fn variant_to_maf(variant: &VariantPosition, decision: &FilterDecision) -> MAFRecord {
    // Select canonical transcript
    let transcript = select_canonical_transcript(&variant.transcripts);

    // Extract gene symbol (use hgnc field which contains the gene symbol)
    let hugo_symbol = transcript
        .as_ref()
        .and_then(|t| t.hgnc.as_deref())
        .unwrap_or("")
        .to_string();

    // Map variant classification
    let variant_classification = transcript
        .as_ref()
        .map(|t| map_variant_classification(&t.consequence))
        .unwrap_or_else(|| "".to_string());

    // Map variant type
    let variant_type = map_variant_type(&variant.variant_type);

    // Extract HGVS notation
    let (hgvsc, hgvsp, hgvsp_short) = extract_hgvs_notation(transcript.as_ref());

    // Transcript ID
    let transcript_id = transcript
        .as_ref()
        .and_then(|t| t.id.as_deref())
        .unwrap_or("")
        .to_string();

    // Extract new annotation fields
    let exon = transcript
        .as_ref()
        .and_then(|t| t.exons.as_deref())
        .unwrap_or("")
        .to_string();

    let consequence = transcript
        .as_ref()
        .map(|t| t.consequence.join(","))
        .unwrap_or_else(|| "".to_string());

    let impact = transcript
        .as_ref()
        .and_then(|t| t.impact.as_deref())
        .map(|s| s.to_uppercase())  // Convert to uppercase (LOW -> LOW, moderate -> MODERATE)
        .unwrap_or_else(|| "".to_string());

    let codons = transcript
        .as_ref()
        .and_then(|t| t.codons.as_deref())
        .unwrap_or("")
        .to_string();

    let amino_acids = transcript
        .as_ref()
        .and_then(|t| t.amino_acids.as_deref())
        .unwrap_or("")
        .to_string();

    let cdna_position = transcript
        .as_ref()
        .and_then(|t| t.cdna_pos.as_deref())
        .unwrap_or("")
        .to_string();

    let cds_position = transcript
        .as_ref()
        .and_then(|t| t.cds_pos.as_deref())
        .unwrap_or("")
        .to_string();

    let protein_position = transcript
        .as_ref()
        .and_then(|t| t.protein_pos.as_deref())
        .unwrap_or("")
        .to_string();

    // dbSNP ID
    let dbsnp_rs = variant
        .dbsnp_ids
        .first()
        .map(|s| s.to_string())
        .unwrap_or_else(|| "".to_string());

    // COSMIC ID
    let cosmic_id = variant
        .cosmic
        .first()
        .and_then(|c| c.id.as_deref())
        .unwrap_or("")
        .to_string();

    // ClinVar information
    let (clinvar_id, clinvar_review_status, clinvar_significance, clinvar_disease) =
        if let Some(entry) = decision
            .primary_evidence
            .as_str()
            .eq("ClinVar")
            .then(|| variant.clinvar.first())
            .flatten()
        {
            (
                entry.id.as_deref().unwrap_or("").to_string(),
                entry.review_status.as_deref().unwrap_or("").to_string(),
                entry.clinical_significance.join(", "),
                entry.phenotypes.join("; "),
            )
        } else {
            (String::new(), String::new(), String::new(), String::new())
        };

    // Predictive scores
    let primate_ai_score = variant
        .primate_ai_3d
        .or(variant.primate_ai)
        .map(|s| format!("{:.4}", s))
        .unwrap_or_else(|| "".to_string());

    let dann_score = variant
        .dann_score
        .map(|s| format!("{:.4}", s))
        .unwrap_or_else(|| "".to_string());

    let revel_score = variant
        .revel_score
        .map(|s| format!("{:.4}", s))
        .unwrap_or_else(|| "".to_string());

    // Population frequency
    let (gnomad_af, gnomad_eas_af) = extract_population_frequencies(variant);

    // Sequencing quality
    let depth = variant
        .total_depth
        .map(|d| d.to_string())
        .unwrap_or_else(|| "".to_string());

    let vaf = variant
        .variant_frequencies
        .as_ref()
        .and_then(|vf| vf.first())
        .map(|v| format!("{:.4}", v))
        .unwrap_or_else(|| "".to_string());

    MAFRecord {
        hugo_symbol,
        chromosome: variant.chromosome.clone(),
        start_position: variant.start,
        end_position: variant.end_pos,
        strand: "+".to_string(),
        variant_classification,
        variant_type,
        reference_allele: variant.reference_allele.clone(),
        tumor_seq_allele1: variant.reference_allele.clone(),
        tumor_seq_allele2: variant.alternate_allele.clone(),
        tumor_sample_barcode: "".to_string(),
        hgvsc,
        hgvsp,
        hgvsp_short,
        transcript_id,
        exon,
        consequence,
        impact,
        codons,
        amino_acids,
        cdna_position,
        cds_position,
        protein_position,
        dbsnp_rs,
        dbsnp_val_status: "".to_string(),
        cosmic_id,
        clinvar_id,
        clinvar_review_status,
        clinvar_significance,
        clinvar_disease,
        primate_ai_score,
        dann_score,
        revel_score,
        gnomad_af,
        gnomad_eas_af,
        depth,
        vaf,
    }
}

fn select_canonical_transcript(
    transcripts: &[TranscriptAnnotation],
) -> Option<TranscriptAnnotation> {
    // Prefer MANE Select transcript
    if let Some(mane) = transcripts
        .iter()
        .find(|t| t.is_mane_select == Some(true))
    {
        return Some(mane.clone());
    }

    // Otherwise, return first transcript
    transcripts.first().cloned()
}

fn map_variant_classification(consequences: &[String]) -> String {
    // Map SO terms to MAF variant classification
    for consequence in consequences {
        let consequence_lower = consequence.to_lowercase();
        let classification = match consequence_lower.as_str() {
            s if s.contains("missense") => "Missense_Mutation",
            s if s.contains("nonsense") || s.contains("stop_gained") => "Nonsense_Mutation",
            s if s.contains("frameshift") => "Frame_Shift_Del",
            s if s.contains("splice_acceptor") || s.contains("splice_donor") => "Splice_Site",
            s if s.contains("inframe_deletion") => "In_Frame_Del",
            s if s.contains("inframe_insertion") => "In_Frame_Ins",
            s if s.contains("start_lost") => "Translation_Start_Site",
            s if s.contains("stop_lost") => "Nonstop_Mutation",
            s if s.contains("synonymous") => "Silent",
            s if s.contains("5_prime_utr") => "5'UTR",
            s if s.contains("3_prime_utr") => "3'UTR",
            s if s.contains("intron") => "Intron",
            _ => continue,
        };
        return classification.to_string();
    }

    "".to_string()
}

fn map_variant_type(variant_type: &str) -> String {
    match variant_type {
        "SNV" => "SNP".to_string(),
        "insertion" => "INS".to_string(),
        "deletion" => "DEL".to_string(),
        "MNV" => "DNP".to_string(),
        _ => variant_type.to_string(),
    }
}

fn extract_hgvs_notation(transcript: Option<&TranscriptAnnotation>) -> (String, String, String) {
    if let Some(t) = transcript {
        let hgvsc = t.hgvsc.as_deref().unwrap_or("").to_string();
        let hgvsp = t.hgvsp.as_deref().unwrap_or("").to_string();

        // Generate short HGVS protein notation (e.g., p.V600E)
        let hgvsp_short = if let Some(hgvsp_str) = &t.hgvsp {
            // Extract short form from full HGVS (e.g., "NP_004324.2:p.Val600Glu" -> "p.V600E")
            if let Some(pos) = hgvsp_str.find(":p.") {
                shorten_hgvsp(&hgvsp_str[pos + 1..])
            } else if hgvsp_str.starts_with("p.") {
                shorten_hgvsp(hgvsp_str)
            } else {
                hgvsp_str.to_string()
            }
        } else {
            String::new()
        };

        (hgvsc, hgvsp, hgvsp_short)
    } else {
        (String::new(), String::new(), String::new())
    }
}

fn shorten_hgvsp(hgvsp: &str) -> String {
    // Convert three-letter amino acid codes to single-letter codes
    hgvsp
        .replace("Ala", "A")
        .replace("Arg", "R")
        .replace("Asn", "N")
        .replace("Asp", "D")
        .replace("Cys", "C")
        .replace("Gln", "Q")
        .replace("Glu", "E")
        .replace("Gly", "G")
        .replace("His", "H")
        .replace("Ile", "I")
        .replace("Leu", "L")
        .replace("Lys", "K")
        .replace("Met", "M")
        .replace("Phe", "F")
        .replace("Pro", "P")
        .replace("Ser", "S")
        .replace("Thr", "T")
        .replace("Trp", "W")
        .replace("Tyr", "Y")
        .replace("Val", "V")
        .replace("Ter", "*")
}

fn extract_population_frequencies(variant: &VariantPosition) -> (String, String) {
    let gnomad_exome = variant
        .population_frequencies
        .iter()
        .find(|pf| pf.source == "gnomad-exome");

    let gnomad_af = gnomad_exome
        .and_then(|pf| pf.all_af)
        .map(|af| format!("{:.6}", af))
        .unwrap_or_else(|| "".to_string());

    let gnomad_eas_af = gnomad_exome
        .and_then(|pf| pf.eas_af)
        .map(|af| format!("{:.6}", af))
        .unwrap_or_else(|| "".to_string());

    (gnomad_af, gnomad_eas_af)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_map_variant_classification() {
        assert_eq!(
            map_variant_classification(&vec!["missense_variant".to_string()]),
            "Missense_Mutation"
        );
        assert_eq!(
            map_variant_classification(&vec!["stop_gained".to_string()]),
            "Nonsense_Mutation"
        );
        assert_eq!(
            map_variant_classification(&vec!["frameshift_variant".to_string()]),
            "Frame_Shift_Del"
        );
        assert_eq!(
            map_variant_classification(&vec!["splice_donor_variant".to_string()]),
            "Splice_Site"
        );
        assert_eq!(
            map_variant_classification(&vec!["synonymous_variant".to_string()]),
            "Silent"
        );
        assert_eq!(
            map_variant_classification(&vec!["inframe_deletion".to_string()]),
            "In_Frame_Del"
        );
    }

    #[test]
    fn test_map_variant_type() {
        assert_eq!(map_variant_type("SNV"), "SNP");
        assert_eq!(map_variant_type("insertion"), "INS");
        assert_eq!(map_variant_type("deletion"), "DEL");
        assert_eq!(map_variant_type("MNV"), "DNP");
    }

    #[test]
    fn test_shorten_hgvsp() {
        assert_eq!(shorten_hgvsp("p.Val600Glu"), "p.V600E");
        assert_eq!(shorten_hgvsp("p.Arg132His"), "p.R132H");
        assert_eq!(shorten_hgvsp("p.Gly12Asp"), "p.G12D");
        assert_eq!(shorten_hgvsp("p.Leu858Arg"), "p.L858R");
    }

    #[test]
    fn test_extract_new_annotation_fields() {
        // Create a test transcript with all new fields
        let transcript = TranscriptAnnotation {
            id: Some("NM_004333.4".to_string()),
            source: Some("RefSeq".to_string()),
            hgnc: Some("BRAF".to_string()),
            consequence: vec!["missense_variant".to_string(), "splice_region_variant".to_string()],
            impact: Some("moderate".to_string()),
            amino_acids: Some("V/E".to_string()),
            cdna_pos: Some("1799/2301".to_string()),
            cds_pos: Some("1799/2301".to_string()),
            exons: Some("15/18".to_string()),
            codons: Some("Gtg/Gag".to_string()),
            protein_pos: Some("600/766".to_string()),
            hgvsc: Some("c.1799T>A".to_string()),
            hgvsp: Some("p.Val600Glu".to_string()),
            is_canonical: Some(true),
            is_mane_select: Some(true),
        };

        // Create a minimal variant for testing
        let variant = VariantPosition {
            chromosome: "chr7".to_string(),
            start: 140453136,
            end_pos: 140453136,
            reference_allele: "A".to_string(),
            alternate_allele: "T".to_string(),
            variant_type: "SNV".to_string(),
            filters: vec!["PASS".to_string()],
            total_depth: Some(100),
            variant_frequencies: Some(vec![0.45]),
            transcripts: vec![transcript],
            clinvar: vec![],
            cosmic: vec![],
            population_frequencies: vec![],
            primate_ai_3d: None,
            primate_ai: None,
            dann_score: None,
            revel_score: None,
            dbsnp_ids: vec![],
        };

        let decision = FilterDecision {
            should_include: true,
            pathogenicity_class: "Pathogenic".to_string(),
            primary_evidence: "ClinVar".to_string(),
            justification: "Test".to_string(),
        };

        let maf = variant_to_maf(&variant, &decision);

        // Verify all new fields are correctly extracted
        assert_eq!(maf.impact, "MODERATE"); // Should be uppercase
        assert_eq!(maf.consequence, "missense_variant,splice_region_variant");
        assert_eq!(maf.exon, "15/18");
        assert_eq!(maf.codons, "Gtg/Gag");
        assert_eq!(maf.amino_acids, "V/E");
        assert_eq!(maf.cdna_position, "1799/2301");
        assert_eq!(maf.cds_position, "1799/2301");
        assert_eq!(maf.protein_position, "600/766");
    }

    #[test]
    fn test_missing_new_fields() {
        // Create a transcript with missing optional fields
        let transcript = TranscriptAnnotation {
            id: Some("NM_001234.1".to_string()),
            source: Some("RefSeq".to_string()),
            hgnc: Some("GENE1".to_string()),
            consequence: vec!["synonymous_variant".to_string()],
            impact: None,
            amino_acids: None,
            cdna_pos: None,
            cds_pos: None,
            exons: None,
            codons: None,
            protein_pos: None,
            hgvsc: None,
            hgvsp: None,
            is_canonical: None,
            is_mane_select: None,
        };

        let variant = VariantPosition {
            chromosome: "chr1".to_string(),
            start: 12345,
            end_pos: 12345,
            reference_allele: "A".to_string(),
            alternate_allele: "G".to_string(),
            variant_type: "SNV".to_string(),
            filters: vec!["PASS".to_string()],
            total_depth: Some(50),
            variant_frequencies: Some(vec![0.3]),
            transcripts: vec![transcript],
            clinvar: vec![],
            cosmic: vec![],
            population_frequencies: vec![],
            primate_ai_3d: None,
            primate_ai: None,
            dann_score: None,
            revel_score: None,
            dbsnp_ids: vec![],
        };

        let decision = FilterDecision {
            should_include: false,
            pathogenicity_class: "Excluded".to_string(),
            primary_evidence: "None".to_string(),
            justification: "Test".to_string(),
        };

        let maf = variant_to_maf(&variant, &decision);

        // Verify missing fields result in empty strings
        assert_eq!(maf.impact, "");
        assert_eq!(maf.exon, "");
        assert_eq!(maf.codons, "");
        assert_eq!(maf.amino_acids, "");
        assert_eq!(maf.cdna_position, "");
        assert_eq!(maf.cds_position, "");
        assert_eq!(maf.protein_position, "");
    }
}
