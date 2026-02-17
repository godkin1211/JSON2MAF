use crate::types::*;
use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use serde_json::Value;
use std::fs::File;
use std::io::{BufReader, Write};
use tempfile::NamedTempFile;

pub fn parse_nirvana_json(file_path: &str) -> Result<(NirvanaHeader, Vec<VariantPosition>)> {
    // First, decompress the file to a temporary location
    let temp_file = NamedTempFile::new().context("Failed to create temporary file")?;
    let temp_json_path = temp_file.path();

    log::info!("Decompressing {} to temporary file...", file_path);
    let file = File::open(file_path)
        .with_context(|| format!("Failed to open input file: {}", file_path))?;

    let reader = BufReader::new(file);
    let mut decoder = MultiGzDecoder::new(reader);

    // Write decompressed data to temporary file
    let bytes_written = std::io::copy(&mut decoder, &mut temp_file.as_file())
        .context("Failed to decompress file")?;

    temp_file.as_file().flush().context("Failed to flush temp file")?;

    log::info!("Decompressed {} bytes", bytes_written);

    log::info!("Decompression complete. Parsing JSON...");

    // Now parse the decompressed JSON file
    let json_file = File::open(temp_json_path)
        .context("Failed to open temporary JSON file")?;
    let json_reader = BufReader::with_capacity(16 * 1024 * 1024, json_file); // 16MB buffer

    // Parse JSON from the decompressed file
    let json: Value = serde_json::from_reader(json_reader)
        .context("Failed to parse JSON from decompressed file")?;

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
        .with_context(|| {
            format!(
                "Failed to parse position. JSON preview: {}",
                serde_json::to_string_pretty(&pos_json).unwrap_or_else(|_| "unable to serialize".to_string()).chars().take(500).collect::<String>()
            )
        })?;

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

    // Extract predictive scores from complex structures
    let primate_ai_3d = variant
        .primate_ai_3d
        .first()
        .and_then(|entry| entry.score);

    let primate_ai = variant
        .primate_ai
        .first()
        .and_then(|entry| entry.score_percentile);

    let dann_score = variant.dann_score;

    let revel_score = variant
        .revel_score
        .as_ref()
        .and_then(|rs| rs.score);

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
        primate_ai_3d,
        primate_ai,
        dann_score,
        revel_score,
        dbsnp_ids: variant.dbsnp.clone(),
    };

    Ok(Some(variant_pos))
}

fn extract_population_frequencies(pos_json: &Value) -> Vec<PopulationFrequency> {
    let mut result = Vec::new();

    // Get first variant from the variants array
    let variant = match pos_json
        .get("variants")
        .and_then(|v| v.as_array())
        .and_then(|arr| arr.first())
    {
        Some(v) => v,
        None => return result,
    };

    // Extract gnomad
    if let Some(gnomad) = variant.get("gnomad") {
        result.push(PopulationFrequency {
            source: "gnomad".to_string(),
            all_af: gnomad.get("allAf").and_then(|v| v.as_f64()),
            eas_af: gnomad.get("easAf").and_then(|v| v.as_f64()),
            afr_af: gnomad.get("afrAf").and_then(|v| v.as_f64()),
            amr_af: gnomad.get("amrAf").and_then(|v| v.as_f64()),
            eur_af: None,
        });
    }

    // Extract gnomad-exome
    if let Some(gnomad_exome) = variant.get("gnomad-exome") {
        result.push(PopulationFrequency {
            source: "gnomad-exome".to_string(),
            all_af: gnomad_exome.get("allAf").and_then(|v| v.as_f64()),
            eas_af: gnomad_exome.get("easAf").and_then(|v| v.as_f64()),
            afr_af: gnomad_exome.get("afrAf").and_then(|v| v.as_f64()),
            amr_af: gnomad_exome.get("amrAf").and_then(|v| v.as_f64()),
            eur_af: None,
        });
    }

    // Extract oneKg (1000 Genomes)
    if let Some(onekg) = variant.get("oneKg") {
        result.push(PopulationFrequency {
            source: "oneKg".to_string(),
            all_af: onekg.get("allAf").and_then(|v| v.as_f64()),
            eas_af: onekg.get("easAf").and_then(|v| v.as_f64()),
            afr_af: onekg.get("afrAf").and_then(|v| v.as_f64()),
            amr_af: onekg.get("amrAf").and_then(|v| v.as_f64()),
            eur_af: onekg.get("eurAf").and_then(|v| v.as_f64()),
        });
    }

    result
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
