use crate::types::*;
use std::collections::HashMap;

pub fn assess_predictive_scores(
    variant: &VariantPosition,
    config: &FilterConfig,
) -> PredictiveAssessment {
    let mut contributing = HashMap::new();
    let mut support_count = 0;

    // Check PrimateAI-3D
    if let Some(primate_ai_3d) = get_primate_ai_score(variant) {
        if primate_ai_3d >= config.min_primate_ai_score {
            contributing.insert("PrimateAI-3D".to_string(), primate_ai_3d);
            support_count += 1;
        }
    }

    // Check REVEL
    if let Some(revel) = get_revel_score(variant) {
        if revel >= config.min_revel_score {
            contributing.insert("REVEL".to_string(), revel);
            support_count += 1;
        }
    }

    // Check DANN
    if let Some(dann) = get_dann_score(variant) {
        if dann >= config.min_dann_score {
            contributing.insert("DANN".to_string(), dann);
            support_count += 1;
        }
    }

    // Check COSMIC (presence indicates positive evidence)
    if is_in_cosmic(variant) {
        contributing.insert("COSMIC".to_string(), 1.0);
        support_count += 1;
    }

    // Determine if should be suggested as likely pathogenic
    let has_primate_ai_3d = contributing.contains_key("PrimateAI-3D");
    let suggests_pathogenic = has_primate_ai_3d || support_count >= 2;

    // Calculate confidence score
    let confidence = calculate_confidence(&contributing, has_primate_ai_3d, support_count);

    PredictiveAssessment {
        suggests_pathogenic,
        contributing_scores: contributing,
        confidence,
        support_count,
        has_primate_ai_support: has_primate_ai_3d,
    }
}

fn calculate_confidence(
    _contributing: &HashMap<String, f64>,
    has_primate_ai: bool,
    support_count: usize,
) -> f64 {
    if support_count == 0 {
        return 0.0;
    }

    // If PrimateAI-3D present, base confidence is 0.7
    let base = if has_primate_ai { 0.7 } else { 0.5 };

    // Each additional supporting score increases confidence
    let confidence = base + (support_count.saturating_sub(1) as f64) * 0.1;

    // Maximum 1.0
    confidence.min(1.0)
}

pub fn get_primate_ai_score(variant: &VariantPosition) -> Option<f64> {
    // Prioritize PrimateAI-3D
    variant.primate_ai_3d.or(variant.primate_ai)
}

pub fn get_dann_score(variant: &VariantPosition) -> Option<f64> {
    variant.dann_score
}

pub fn get_revel_score(variant: &VariantPosition) -> Option<f64> {
    variant.revel_score
}

pub fn is_in_cosmic(variant: &VariantPosition) -> bool {
    !variant.cosmic.is_empty()
}

pub fn count_supporting_predictive_scores(assessment: &PredictiveAssessment) -> usize {
    assessment.support_count
}

pub fn has_primate_ai_support(assessment: &PredictiveAssessment) -> bool {
    assessment.has_primate_ai_support
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_primate_ai_alone_supports() {
        let mut variant = create_test_variant();
        variant.primate_ai_3d = Some(0.85);

        let config = FilterConfig::default();
        let assessment = assess_predictive_scores(&variant, &config);

        assert!(assessment.suggests_pathogenic);
        assert!(assessment.has_primate_ai_support);
        assert_eq!(assessment.support_count, 1);
    }

    #[test]
    fn test_two_scores_support() {
        let mut variant = create_test_variant();
        variant.revel_score = Some(0.8);
        variant.dann_score = Some(0.97);

        let config = FilterConfig::default();
        let assessment = assess_predictive_scores(&variant, &config);

        assert!(assessment.suggests_pathogenic);
        assert!(!assessment.has_primate_ai_support);
        assert_eq!(assessment.support_count, 2);
    }

    #[test]
    fn test_insufficient_support() {
        let mut variant = create_test_variant();
        variant.revel_score = Some(0.8);

        let config = FilterConfig::default();
        let assessment = assess_predictive_scores(&variant, &config);

        assert!(!assessment.suggests_pathogenic);
        assert_eq!(assessment.support_count, 1);
    }

    fn create_test_variant() -> VariantPosition {
        VariantPosition {
            chromosome: "chr1".to_string(),
            start: 100,
            end_pos: 100,
            reference_allele: "A".to_string(),
            alternate_allele: "T".to_string(),
            variant_type: "SNV".to_string(),
            filters: vec!["PASS".to_string()],
            total_depth: Some(50),
            variant_frequencies: Some(vec![0.05]),
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
