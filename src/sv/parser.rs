/// SV JSON parser using the same per-position `serde_json::Value` extraction
/// approach as before, but driven by a single-pass stream instead of a
/// whole-file DOM. Each `positions[]` element is deserialized as its own
/// small, bounded `Value` and dropped immediately after processing, so
/// memory stays proportional to one SV record rather than the whole
/// (potentially multi-GB decompressed) file.
///
/// Key differences from SNV parsing:
/// - Reads `svEnd` and `svLength` from position level (not inside variants)
/// - Reads `clingen` array from position level
/// - Sample read support from `splitReadCounts` and `pairedEndReadCounts` (integer arrays)
/// - Classifies SV type from altAllele string; skips BND (breakend) variants
use anyhow::{Context, Result};
use serde_json::Value;

use crate::json_stream::stream_positions;
use crate::{NirvanaHeader, TranscriptAnnotation};
use super::types::*;

// ============================================================================
// Public entry point
// ============================================================================

/// Streams `SVPosition` records straight out of the gzipped Nirvana SV JSON,
/// one at a time. `on_position` is invoked once per surviving position
/// (BND / unsupported-type positions are skipped, and per-position parse
/// errors are logged and skipped, matching the previous behavior).
pub fn parse_sv_nirvana_streaming<F>(file_path: &str, mut on_position: F) -> Result<NirvanaHeader>
where
    F: FnMut(SVPosition) -> Result<()>,
{
    // `header` and `positions` callbacks are both live for the duration of the
    // stream, so `sample_name` needs interior mutability to be written by one
    // and read by the other.
    let sample_name = std::cell::RefCell::new(String::new());

    stream_positions::<Value, NirvanaHeader, _, _>(
        file_path,
        |header| {
            *sample_name.borrow_mut() = header.samples.first().cloned().unwrap_or_default();
            Ok(())
        },
        |pos_json: Value| {
            match parse_sv_position(&pos_json, &sample_name.borrow()) {
                Ok(Some(sv_pos)) => on_position(sv_pos)?,
                Ok(None) => {} // Skipped (BND or unsupported type)
                Err(e) => log::warn!("Skipping position due to parse error: {}", e),
            }
            Ok(())
        },
    )
    .context("Failed to parse SV JSON")
}

/// Convenience wrapper that collects every SV position into a `Vec`. Kept
/// for callers that want the whole result at once; the CLI binary
/// (`json2sv`) uses `parse_sv_nirvana_streaming` directly so it never has
/// to hold every SV position in memory at once.
pub fn parse_sv_nirvana_json(file_path: &str) -> Result<(NirvanaHeader, Vec<SVPosition>)> {
    let mut sv_positions = Vec::new();
    let header = parse_sv_nirvana_streaming(file_path, |sv_pos| {
        sv_positions.push(sv_pos);
        Ok(())
    })?;
    Ok((header, sv_positions))
}

// ============================================================================
// Position parsing
// ============================================================================

fn parse_sv_position(pos: &Value, sample_name: &str) -> Result<Option<SVPosition>> {
    let chromosome = pos["chromosome"]
        .as_str()
        .ok_or_else(|| anyhow::anyhow!("Missing chromosome"))?
        .to_string();

    let start = pos["position"]
        .as_i64()
        .ok_or_else(|| anyhow::anyhow!("Missing position at {}", chromosome))?;

    let alt_alleles = pos["altAlleles"]
        .as_array()
        .ok_or_else(|| anyhow::anyhow!("Missing altAlleles at {}:{}", chromosome, start))?;

    // Only process the first alt allele; skip if BND
    let alt = alt_alleles
        .first()
        .and_then(|v| v.as_str())
        .ok_or_else(|| anyhow::anyhow!("Empty altAlleles at {}:{}", chromosome, start))?;

    let ref_allele = pos["refAllele"]
        .as_str()
        .unwrap_or("N")
        .to_string();

    // Classify SV type; returns None for BND or same-length alleles
    let sv_classification = match classify_sv(alt, &ref_allele, pos) {
        Some(c) => c,
        None => return Ok(None),
    };

    let filters = parse_string_array(&pos["filters"]);
    let read_support = parse_read_support(pos);
    let transcripts = parse_sv_transcripts(pos);
    let clingen = parse_clingen(pos);

    Ok(Some(SVPosition {
        chromosome,
        start,
        end_pos: sv_classification.end_pos,
        sv_type: sv_classification.sv_type,
        sv_length: sv_classification.sv_length,
        is_symbolic: sv_classification.is_symbolic,
        reference_allele: ref_allele,
        alternate_allele: alt.to_string(),
        filters,
        read_support,
        transcripts,
        clingen,
        sample_name: sample_name.to_string(),
    }))
}

// ============================================================================
// SV classification logic
// ============================================================================

struct SVClassification {
    sv_type: SVType,
    is_symbolic: bool,
    sv_length: i64,
    end_pos: i64,
}

/// Returns None for BND variants (alt contains '[' or ']') or equal-length alleles.
fn classify_sv(alt: &str, reference: &str, pos: &Value) -> Option<SVClassification> {
    // Skip breakend variants
    if alt.contains('[') || alt.contains(']') {
        return None;
    }

    match alt {
        "<DEL>" => {
            let end_pos = pos["svEnd"].as_i64().unwrap_or_else(|| {
                pos["position"].as_i64().unwrap_or(0)
            });
            let sv_length = pos["svLength"].as_i64().unwrap_or(end_pos - pos["position"].as_i64().unwrap_or(0)).abs();
            Some(SVClassification {
                sv_type: SVType::Del,
                is_symbolic: true,
                sv_length,
                end_pos,
            })
        }
        "<INS>" => {
            let end_pos = pos["svEnd"].as_i64().unwrap_or_else(|| {
                pos["position"].as_i64().unwrap_or(0)
            });
            let sv_length = pos["svLength"].as_i64().unwrap_or(0).abs();
            Some(SVClassification {
                sv_type: SVType::Ins,
                is_symbolic: true,
                sv_length,
                end_pos,
            })
        }
        _ => {
            // Sequence-resolved: infer type from allele lengths
            let ref_len = reference.len() as i64;
            let alt_len = alt.len() as i64;

            if ref_len > alt_len {
                // Deletion: reference is longer
                let sv_length = ref_len - alt_len;
                let end_pos = pos["position"].as_i64().unwrap_or(0) + ref_len - 1;
                Some(SVClassification {
                    sv_type: SVType::Del,
                    is_symbolic: false,
                    sv_length,
                    end_pos,
                })
            } else if alt_len > ref_len {
                // Insertion: alt is longer
                let sv_length = alt_len - ref_len;
                let end_pos = pos["position"].as_i64().unwrap_or(0);
                Some(SVClassification {
                    sv_type: SVType::Ins,
                    is_symbolic: false,
                    sv_length,
                    end_pos,
                })
            } else {
                // Same length — not a DEL/INS (could be MNV), skip
                None
            }
        }
    }
}

// ============================================================================
// Read support extraction
// ============================================================================

/// Extract read counts from `samples[0]`.
/// SV samples use `splitReadCounts` and `pairedEndReadCounts` (integer arrays),
/// unlike SNV samples which use `totalDepth` and `variantFrequencies` (floats).
fn parse_read_support(pos: &Value) -> SVReadSupport {
    let sample = match pos["samples"].as_array().and_then(|a| a.first()) {
        Some(s) => s,
        None => return SVReadSupport::default(),
    };

    let split_read_alt = sample["splitReadCounts"][0].as_i64().unwrap_or(0);
    let split_read_ref = sample["splitReadCounts"][1].as_i64().unwrap_or(0);
    let paired_end_alt = sample["pairedEndReadCounts"][0].as_i64().unwrap_or(0);
    let paired_end_ref = sample["pairedEndReadCounts"][1].as_i64().unwrap_or(0);

    SVReadSupport {
        split_read_alt,
        split_read_ref,
        paired_end_alt,
        paired_end_ref,
    }
}

// ============================================================================
// Transcript parsing (reuses TranscriptAnnotation from types.rs)
// ============================================================================

/// Extracts transcript annotations from `variants[0].transcripts`.
/// Reuses `serde_json::from_value::<TranscriptAnnotation>` — the same struct
/// used for SNVs works here because all its fields are optional.
fn parse_sv_transcripts(pos: &Value) -> Vec<TranscriptAnnotation> {
    let transcripts_json = match pos["variants"][0]["transcripts"].as_array() {
        Some(a) => a,
        None => return Vec::new(),
    };

    transcripts_json
        .iter()
        .filter_map(|t| serde_json::from_value::<TranscriptAnnotation>(t.clone()).ok())
        .collect()
}

// ============================================================================
// ClinGen parsing
// ============================================================================

fn parse_clingen(pos: &Value) -> Vec<ClinGenEntry> {
    let clingen_json = match pos["clingen"].as_array() {
        Some(a) => a,
        None => return Vec::new(),
    };

    clingen_json
        .iter()
        .map(|c| ClinGenEntry {
            id: c["id"].as_str().unwrap_or("").to_string(),
            variant_type: c["variantType"].as_str().unwrap_or("").to_string(),
            clinical_interpretation: c["clinicalInterpretation"]
                .as_str()
                .unwrap_or("")
                .to_string(),
            phenotypes: parse_string_array(&c["phenotypes"]),
        })
        .collect()
}

// ============================================================================
// Helpers
// ============================================================================

fn parse_string_array(value: &Value) -> Vec<String> {
    value
        .as_array()
        .map(|arr| {
            arr.iter()
                .filter_map(|v| v.as_str().map(str::to_string))
                .collect()
        })
        .unwrap_or_default()
}
