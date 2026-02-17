use crate::types::*;

pub fn apply_quality_filters(
    variant: &VariantPosition,
    config: &FilterConfig,
) -> QualityFilterResult {
    // Check VCF filters field: only accept ["PASS"]
    if !(variant.filters.len() == 1 && variant.filters[0] == "PASS") {
        let filters_str = variant.filters.join(", ");
        return QualityFilterResult {
            passes_quality: false,
            failure_reason: Some(format!("Failed VCF filters: [{}]", filters_str)),
            depth: variant.total_depth,
            variant_frequency: get_variant_frequency(variant),
            eas_allele_frequency: None,
        };
    }

    // Check sequencing quality
    if let Some((false, reason)) = check_sequencing_quality(variant, config) {
        return QualityFilterResult {
            passes_quality: false,
            failure_reason: Some(reason),
            depth: variant.total_depth,
            variant_frequency: get_variant_frequency(variant),
            eas_allele_frequency: None,
        };
    }

    // Check population frequency
    let (pop_pass, pop_reason, eas_af) = check_population_frequency(variant, config);
    if !pop_pass {
        return QualityFilterResult {
            passes_quality: false,
            failure_reason: pop_reason,
            depth: variant.total_depth,
            variant_frequency: get_variant_frequency(variant),
            eas_allele_frequency: eas_af,
        };
    }

    // All passed
    QualityFilterResult {
        passes_quality: true,
        failure_reason: None,
        depth: variant.total_depth,
        variant_frequency: get_variant_frequency(variant),
        eas_allele_frequency: eas_af,
    }
}

fn check_sequencing_quality(
    variant: &VariantPosition,
    config: &FilterConfig,
) -> Option<(bool, String)> {
    // Check sequencing depth
    let depth = variant.total_depth?;
    if depth < config.min_total_depth {
        return Some((
            false,
            format!(
                "Low sequencing depth ({} < {})",
                depth, config.min_total_depth
            ),
        ));
    }

    // Check variant frequency
    let vaf = get_variant_frequency(variant)?;
    if vaf < config.min_variant_frequency {
        return Some((
            false,
            format!(
                "Low variant frequency ({:.4} < {})",
                vaf, config.min_variant_frequency
            ),
        ));
    }

    Some((true, String::new()))
}

fn check_population_frequency(
    variant: &VariantPosition,
    config: &FilterConfig,
) -> (bool, Option<String>, Option<f64>) {
    // Try to extract easAf from gnomad-exome
    if let Some(gnomad_af) = extract_gnomad_exome_eas_af(variant) {
        if gnomad_af > config.max_eas_af {
            return (
                false,
                Some(format!(
                    "High East Asian AF in gnomAD-exome ({:.4} > {})",
                    gnomad_af, config.max_eas_af
                )),
                Some(gnomad_af),
            );
        }
        return (true, None, Some(gnomad_af));
    }

    // Try to extract easAf from oneKg
    if let Some(onekg_af) = extract_onekg_eas_af(variant) {
        if onekg_af > config.max_eas_af {
            return (
                false,
                Some(format!(
                    "High East Asian AF in 1000G ({:.4} > {})",
                    onekg_af, config.max_eas_af
                )),
                Some(onekg_af),
            );
        }
        return (true, None, Some(onekg_af));
    }

    // No population frequency data, consider as pass (conservative strategy)
    (true, None, None)
}

fn extract_gnomad_exome_eas_af(variant: &VariantPosition) -> Option<f64> {
    variant
        .population_frequencies
        .iter()
        .find(|pf| pf.source == "gnomad-exome")
        .and_then(|pf| pf.eas_af)
}

fn extract_onekg_eas_af(variant: &VariantPosition) -> Option<f64> {
    variant
        .population_frequencies
        .iter()
        .find(|pf| pf.source == "oneKg")
        .and_then(|pf| pf.eas_af)
}

fn get_variant_frequency(variant: &VariantPosition) -> Option<f64> {
    variant
        .variant_frequencies
        .as_ref()
        .and_then(|vf| vf.first().copied())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quality_filter_pass() {
        let variant = create_test_variant(50, 0.05);
        let config = FilterConfig::default();
        let result = apply_quality_filters(&variant, &config);
        assert!(result.passes_quality);
    }

    #[test]
    fn test_quality_filter_low_depth() {
        let variant = create_test_variant(20, 0.05);
        let config = FilterConfig::default();
        let result = apply_quality_filters(&variant, &config);
        assert!(!result.passes_quality);
        assert!(result
            .failure_reason
            .unwrap()
            .contains("Low sequencing depth"));
    }

    #[test]
    fn test_quality_filter_low_vaf() {
        let variant = create_test_variant(50, 0.01);
        let config = FilterConfig::default();
        let result = apply_quality_filters(&variant, &config);
        assert!(!result.passes_quality);
        assert!(result
            .failure_reason
            .unwrap()
            .contains("Low variant frequency"));
    }

    fn create_test_variant(depth: i32, vaf: f64) -> VariantPosition {
        VariantPosition {
            chromosome: "chr1".to_string(),
            start: 100,
            end_pos: 100,
            reference_allele: "A".to_string(),
            alternate_allele: "T".to_string(),
            variant_type: "SNV".to_string(),
            filters: vec!["PASS".to_string()],
            total_depth: Some(depth),
            variant_frequencies: Some(vec![vaf]),
            transcripts: vec![],
            clinvar: vec![],
            cosmic: vec![],
            population_frequencies: vec![],
            primate_ai_3d: None,
            primate_ai: None,
            dann_score: None,
            revel_score: None,
            dbsnp_ids: vec![],
        }
    }
}
