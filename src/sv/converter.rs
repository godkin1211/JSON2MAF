/// Converts `SVPosition` to `SVRecord` (output TSV row).
///
/// Key decisions made here:
/// - **Gene symbol**: collected from all transcripts, deduplicated, sorted
/// - **Consequence**: most severe term wins (priority list mirrors VEP severity)
/// - **Transcript**: MANE Select preferred; fallback to first transcript with HGVSc
/// - **ClinGen**: pathogenic entry preferred over first available
use std::collections::BTreeSet;

use crate::TranscriptAnnotation;
use super::types::*;

/// Consequence severity ranking (most severe first).
/// Mirrors the order used by Ensembl VEP for SV consequences.
const CONSEQUENCE_SEVERITY: &[&str] = &[
    "transcript_ablation",
    "frameshift_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "inframe_deletion",
    "inframe_insertion",
    "feature_truncation",
    "feature_elongation",
    "exon_loss_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "transcript_variant",
    "downstream_gene_variant",
    "upstream_gene_variant",
    "intergenic_variant",
];

pub fn sv_position_to_record(pos: &SVPosition) -> SVRecord {
    let hugo_symbol = collect_gene_symbols(&pos.transcripts);
    let best_tx = select_best_transcript(&pos.transcripts);

    let all_consequences: Vec<String> = pos.transcripts
        .iter()
        .flat_map(|t| t.consequence.iter().cloned())
        .collect();
    let variant_classification = pick_most_severe_consequence(&all_consequences);

    let (clingen_id, clingen_interpretation, clingen_phenotypes) =
        select_clingen_entry(&pos.clingen);

    let total_alt = pos.read_support.total_alt();
    let total_ref = pos.read_support.total_ref();
    let vaf = pos.read_support.vaf()
        .map(|v| format!("{:.4}", v))
        .unwrap_or_else(|| ".".to_string());

    SVRecord {
        hugo_symbol,
        chromosome: pos.chromosome.clone(),
        start_position: pos.start,
        end_position: pos.end_pos,
        sv_type: pos.sv_type.as_str().to_string(),
        sv_length: pos.sv_length,
        variant_classification,
        hgvsc: best_tx.and_then(|t| t.hgvsc.clone()).unwrap_or_default(),
        hgvsp: best_tx.and_then(|t| t.hgvsp.clone()).unwrap_or_default(),
        transcript_id: best_tx.and_then(|t| t.id.clone()).unwrap_or_default(),
        split_read_alt: pos.read_support.split_read_alt,
        split_read_ref: pos.read_support.split_read_ref,
        paired_end_alt: pos.read_support.paired_end_alt,
        paired_end_ref: pos.read_support.paired_end_ref,
        total_alt_support: total_alt,
        total_ref_support: total_ref,
        vaf,
        filters: pos.filters.join(";"),
        tumor_sample_barcode: pos.sample_name.clone(),
        clingen_id,
        clingen_interpretation,
        clingen_phenotypes,
    }
}

// ============================================================================
// Helpers
// ============================================================================

/// Collects unique gene symbols from transcripts, sorted alphabetically.
/// Uses BTreeSet to guarantee deterministic order without post-sort.
fn collect_gene_symbols(transcripts: &[TranscriptAnnotation]) -> String {
    let genes: BTreeSet<String> = transcripts
        .iter()
        .filter_map(|t| t.hgnc.clone())
        .filter(|g| !g.is_empty())
        .collect();
    genes.into_iter().collect::<Vec<_>>().join(";")
}

/// Selects the most relevant transcript for HGVSc/HGVSp annotation.
/// Priority: MANE Select with HGVSc > MANE Select > any with HGVSc > first
fn select_best_transcript(transcripts: &[TranscriptAnnotation]) -> Option<&TranscriptAnnotation> {
    transcripts
        .iter()
        .find(|t| t.is_mane_select == Some(true) && t.hgvsc.is_some())
        .or_else(|| transcripts.iter().find(|t| t.is_mane_select == Some(true)))
        .or_else(|| transcripts.iter().find(|t| t.hgvsc.is_some()))
        .or_else(|| transcripts.first())
}

/// Returns the most severe consequence term found in the list.
fn pick_most_severe_consequence(consequences: &[String]) -> String {
    for &severe in CONSEQUENCE_SEVERITY {
        if consequences.iter().any(|c| c == severe) {
            return severe.to_string();
        }
    }
    consequences.first().cloned().unwrap_or_default()
}

/// Selects the most clinically relevant ClinGen entry.
/// Prefers pathogenic over other interpretations.
fn select_clingen_entry(clingen: &[ClinGenEntry]) -> (String, String, String) {
    let selected = clingen
        .iter()
        .find(|c| c.clinical_interpretation == "pathogenic")
        .or_else(|| clingen.first());

    if let Some(c) = selected {
        (
            c.id.clone(),
            c.clinical_interpretation.clone(),
            c.phenotypes.join("|"),
        )
    } else {
        (String::new(), String::new(), String::new())
    }
}
