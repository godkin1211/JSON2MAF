use crate::types::*;
use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use serde_json::Value;
use std::fs::File;
use std::io::{BufReader, Read};

pub fn parse_nirvana_json(file_path: &str) -> Result<(NirvanaHeader, Vec<VariantPosition>)> {
    let file = File::open(file_path)
        .with_context(|| format!("Failed to open input file: {}", file_path))?;

    let reader = BufReader::new(file);
    let mut decoder = GzDecoder::new(reader);
    let mut content = String::new();
    decoder
        .read_to_string(&mut content)
        .context("Failed to decompress gzip file")?;

    let json: Value = serde_json::from_str(&content).context("Failed to parse JSON")?;

    // Parse header
    let header: NirvanaHeader = serde_json::from_value(
        json.get("header")
            .ok_or_else(|| anyhow::anyhow!("No header found in JSON"))?
            .clone(),
    )
    .context("Failed to parse header")?;

    // Parse positions
    let positions_json = json
        .get("positions")
        .ok_or_else(|| anyhow::anyhow!("No positions found in JSON"))?
        .as_array()
        .ok_or_else(|| anyhow::anyhow!("Positions is not an array"))?;

    let mut variant_positions = Vec::new();

    for pos_json in positions_json {
        if let Some(variant_pos) = parse_position(pos_json)? {
            variant_positions.push(variant_pos);
        }
    }

    Ok((header, variant_positions))
}

fn parse_position(pos_json: &Value) -> Result<Option<VariantPosition>> {
    let position: Position = serde_json::from_value(pos_json.clone())
        .context("Failed to parse position")?;

    // Check if there are variants
    if position.variants.is_empty() {
        return Ok(None);
    }

    let variant = &position.variants[0];

    // Get filters or default to ["PASS"]
    let filters = if position.filters.is_empty() {
        vec!["PASS".to_string()]
    } else {
        position.filters.clone()
    };

    // Get sample information
    let (total_depth, variant_frequencies) = if let Some(sample) = position.samples.first() {
        (sample.total_depth, sample.variant_frequencies.clone())
    } else {
        (None, None)
    };

    // Get alternate allele
    let alternate_allele = position
        .alternate_alleles
        .first()
        .ok_or_else(|| anyhow::anyhow!("No alternate alleles found"))?
        .clone();

    // Parse population frequencies (would need to be extracted from a different part of JSON)
    // For now, we'll use an empty vec - this should be populated based on actual JSON structure
    let population_frequencies = Vec::new();

    let variant_pos = VariantPosition {
        chromosome: position.chromosome,
        start: position.position,
        end_pos: position.position,
        reference_allele: position.reference_allele,
        alternate_allele,
        variant_type: variant.variant_type.clone(),
        filters,
        total_depth,
        variant_frequencies,
        transcripts: variant.transcripts.clone(),
        clinvar: variant.clinvar.clone(),
        cosmic: variant.cosmic.clone(),
        population_frequencies,
        primate_ai_3d: variant.primate_ai_3d,
        primate_ai: variant.primate_ai,
        dann_score: variant.dann_score,
        revel_score: variant.revel_score,
        dbsnp_ids: variant.dbsnp.clone(),
    };

    Ok(Some(variant_pos))
}

pub fn parse_nirvana_streaming<F>(
    file_path: &str,
    _config: &FilterConfig,
    mut process_fn: F,
) -> Result<NirvanaHeader>
where
    F: FnMut(&VariantPosition, usize) -> Result<()>,
{
    let (header, variants) = parse_nirvana_json(file_path)?;

    for (idx, variant) in variants.iter().enumerate() {
        process_fn(variant, idx)?;
    }

    Ok(header)
}
