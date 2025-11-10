use crate::types::*;
use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use serde_json::Value;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub fn parse_nirvana_json(file_path: &str) -> Result<(NirvanaHeader, Vec<VariantPosition>)> {
    let file = File::open(file_path)
        .with_context(|| format!("Failed to open input file: {}", file_path))?;

    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);

    let mut lines = reader.lines();

    // Parse header from first line
    let header_line = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("Empty JSON file"))?
        .context("Failed to read header line")?;

    let header_json: Value = serde_json::from_str(&header_line)
        .context("Failed to parse header JSON")?;

    let header: NirvanaHeader = serde_json::from_value(header_json)
        .context("Failed to parse header structure")?;

    // Parse positions from subsequent lines
    let mut variant_positions = Vec::new();

    for (line_num, line_result) in lines.enumerate() {
        let line = line_result
            .with_context(|| format!("Failed to read line {}", line_num + 2))?;

        // Skip empty lines
        if line.trim().is_empty() {
            continue;
        }

        let pos_json: Value = serde_json::from_str(&line)
            .with_context(|| format!("Failed to parse JSON at line {}", line_num + 2))?;

        if let Some(variant_pos) = parse_position(&pos_json)? {
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

    // Extract population frequencies from the variant
    let population_frequencies = extract_population_frequencies(&pos_json);

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

fn extract_population_frequencies(pos_json: &Value) -> Vec<PopulationFrequency> {
    // Try to get population frequencies from the variants array
    if let Some(variants) = pos_json.get("variants").and_then(|v| v.as_array()) {
        if let Some(variant) = variants.first() {
            if let Some(pop_freqs) = variant.get("populationFrequencies").and_then(|p| p.as_array()) {
                return pop_freqs
                    .iter()
                    .filter_map(|pf| {
                        serde_json::from_value::<PopulationFrequency>(pf.clone()).ok()
                    })
                    .collect();
            }
        }
    }
    Vec::new()
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
