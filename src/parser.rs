use crate::json_stream::stream_positions;
use crate::types::*;
use anyhow::Result;

/// Streams `Position` records straight out of the gzipped Nirvana JSON,
/// one at a time, without ever holding the whole file (or even the whole
/// `positions` array) in memory as a `serde_json::Value` DOM. `on_position`
/// is invoked once per position in file order; return `Err` from it to
/// abort parsing early.
pub fn parse_nirvana_streaming<F>(file_path: &str, mut on_position: F) -> Result<NirvanaHeader>
where
    F: FnMut(Position) -> Result<()>,
{
    stream_positions::<Position, NirvanaHeader, _, _>(file_path, |_header| Ok(()), &mut on_position)
}

/// Convenience wrapper that collects every position into a `Vec`. Kept for
/// callers (tests, small files) that want the whole result at once; the
/// real CLI pipeline in `main.rs` uses `parse_nirvana_streaming` directly so
/// it never has to hold the full variant set in memory.
pub fn parse_nirvana_json(file_path: &str) -> Result<(NirvanaHeader, Vec<VariantPosition>)> {
    let mut variant_positions = Vec::new();

    let header = parse_nirvana_streaming(file_path, |position| {
        if let Some(variant_pos) = position_to_variant(position)? {
            variant_positions.push(variant_pos);
        }
        Ok(())
    })?;

    Ok((header, variant_positions))
}

/// Converts one already-deserialized `Position` into a `VariantPosition`,
/// or `None` if the position has no variants (nothing to report).
///
/// Takes `Position` by value so fields can be moved into the result instead
/// of cloned — this used to go through a `serde_json::Value` clone per
/// position, which was pure overhead once `Position` is deserialized
/// directly from the JSON stream.
pub fn position_to_variant(position: Position) -> Result<Option<VariantPosition>> {
    if position.variants.is_empty() {
        return Ok(None);
    }

    let variant = position
        .variants
        .into_iter()
        .next()
        .expect("checked non-empty above");

    let filters = if position.filters.is_empty() {
        vec!["PASS".to_string()]
    } else {
        position.filters
    };

    let (total_depth, variant_frequencies) = match position.samples.into_iter().next() {
        Some(sample) => (sample.total_depth, sample.variant_frequencies),
        None => (None, None),
    };

    let alternate_allele = position
        .alternate_alleles
        .into_iter()
        .next()
        .ok_or_else(|| anyhow::anyhow!("No alternate alleles found"))?;

    let population_frequencies = extract_population_frequencies(&variant);

    let primate_ai_3d = variant.primate_ai_3d.first().and_then(|entry| entry.score);
    let primate_ai = variant
        .primate_ai
        .first()
        .and_then(|entry| entry.score_percentile);
    let dann_score = variant.dann_score;
    let revel_score = variant.revel_score.as_ref().and_then(|rs| rs.score);

    Ok(Some(VariantPosition {
        chromosome: position.chromosome,
        start: position.position,
        end_pos: position.position,
        reference_allele: position.reference_allele,
        alternate_allele,
        variant_type: variant.variant_type,
        filters,
        total_depth,
        variant_frequencies,
        transcripts: variant.transcripts,
        clinvar: variant.clinvar,
        cosmic: variant.cosmic,
        population_frequencies,
        primate_ai_3d,
        primate_ai,
        dann_score,
        revel_score,
        dbsnp_ids: variant.dbsnp,
    }))
}

fn extract_population_frequencies(variant: &Variant) -> Vec<PopulationFrequency> {
    let mut result = Vec::new();

    if let Some(gnomad) = &variant.gnomad {
        result.push(PopulationFrequency {
            source: "gnomad".to_string(),
            all_af: gnomad.all_af,
            eas_af: gnomad.eas_af,
            afr_af: gnomad.afr_af,
            amr_af: gnomad.amr_af,
            eur_af: None,
        });
    }

    if let Some(gnomad_exome) = &variant.gnomad_exome {
        result.push(PopulationFrequency {
            source: "gnomad-exome".to_string(),
            all_af: gnomad_exome.all_af,
            eas_af: gnomad_exome.eas_af,
            afr_af: gnomad_exome.afr_af,
            amr_af: gnomad_exome.amr_af,
            eur_af: None,
        });
    }

    if let Some(onekg) = &variant.one_kg {
        result.push(PopulationFrequency {
            source: "oneKg".to_string(),
            all_af: onekg.all_af,
            eas_af: onekg.eas_af,
            afr_af: onekg.afr_af,
            amr_af: onekg.amr_af,
            eur_af: onekg.eur_af,
        });
    }

    result
}
