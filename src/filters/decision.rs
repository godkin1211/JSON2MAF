use crate::types::*;

pub fn make_filter_decision(
    _variant: &VariantPosition,
    clinvar_assessment: &ClinVarAssessment,
    predictive_assessment: &PredictiveAssessment,
) -> FilterDecision {
    make_filter_decision_with_config(_variant, clinvar_assessment, predictive_assessment, false)
}

pub fn make_filter_decision_with_config(
    _variant: &VariantPosition,
    clinvar_assessment: &ClinVarAssessment,
    predictive_assessment: &PredictiveAssessment,
    exclude_benign: bool,
) -> FilterDecision {
    // Priority 1: ClinVar Pathogenic (takes precedence over benign)
    if clinvar_assessment.is_pathogenic {
        return FilterDecision {
            should_include: true,
            pathogenicity_class: "Pathogenic".to_string(),
            primary_evidence: "ClinVar".to_string(),
            justification: format!(
                "ClinVar pathogenic variant (confidence: {})",
                clinvar_assessment.confidence_level
            ),
        };
    }

    // Priority 2: ClinVar Likely Pathogenic
    if clinvar_assessment.is_likely_pathogenic {
        return FilterDecision {
            should_include: true,
            pathogenicity_class: "Likely pathogenic".to_string(),
            primary_evidence: "ClinVar".to_string(),
            justification: format!(
                "ClinVar likely pathogenic variant (confidence: {})",
                clinvar_assessment.confidence_level
            ),
        };
    }

    // Check for benign variants before considering predictive scores
    // Only filter benign if exclude_benign is enabled AND no pathogenic evidence from ClinVar
    if exclude_benign && (clinvar_assessment.is_benign || clinvar_assessment.is_likely_benign) {
        let benign_class = if clinvar_assessment.is_benign {
            "Benign"
        } else {
            "Likely benign"
        };
        return FilterDecision {
            should_include: false,
            pathogenicity_class: "Excluded (Benign)".to_string(),
            primary_evidence: "ClinVar".to_string(),
            justification: format!(
                "ClinVar {} variant (confidence: {})",
                benign_class, clinvar_assessment.confidence_level
            ),
        };
    }

    // Priority 3: Predictive scores suggest pathogenic
    if predictive_assessment.suggests_pathogenic {
        let score_names: Vec<String> = predictive_assessment
            .contributing_scores
            .keys()
            .cloned()
            .collect();

        return FilterDecision {
            should_include: true,
            pathogenicity_class: "Likely pathogenic".to_string(),
            primary_evidence: "Predictive".to_string(),
            justification: format!(
                "Supported by predictive scores: {} (confidence: {:.2})",
                score_names.join(", "),
                predictive_assessment.confidence
            ),
        };
    }

    // Exclude variant
    FilterDecision {
        should_include: false,
        pathogenicity_class: "Excluded".to_string(),
        primary_evidence: "None".to_string(),
        justification: "Insufficient evidence for pathogenicity".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_clinvar_pathogenic_decision() {
        let variant = create_test_variant();
        let clinvar = ClinVarAssessment {
            is_pathogenic: true,
            is_likely_pathogenic: false,
            is_benign: false,
            is_likely_benign: false,
            selected_entry: None,
            confidence_level: "high".to_string(),
            reason: "ClinVar pathogenic".to_string(),
        };
        let predictive = create_empty_predictive();

        let decision = make_filter_decision(&variant, &clinvar, &predictive);

        assert!(decision.should_include);
        assert_eq!(decision.pathogenicity_class, "Pathogenic");
        assert_eq!(decision.primary_evidence, "ClinVar");
    }

    #[test]
    fn test_predictive_scores_decision() {
        let variant = create_test_variant();
        let clinvar = create_empty_clinvar();
        let predictive = PredictiveAssessment {
            suggests_pathogenic: true,
            contributing_scores: {
                let mut scores = HashMap::new();
                scores.insert("REVEL".to_string(), 0.8);
                scores.insert("DANN".to_string(), 0.97);
                scores
            },
            confidence: 0.7,
            support_count: 2,
            has_primate_ai_support: false,
        };

        let decision = make_filter_decision(&variant, &clinvar, &predictive);

        assert!(decision.should_include);
        assert_eq!(decision.pathogenicity_class, "Likely pathogenic");
        assert_eq!(decision.primary_evidence, "Predictive");
    }

    #[test]
    fn test_exclude_decision() {
        let variant = create_test_variant();
        let clinvar = create_empty_clinvar();
        let predictive = create_empty_predictive();

        let decision = make_filter_decision(&variant, &clinvar, &predictive);

        assert!(!decision.should_include);
        assert_eq!(decision.pathogenicity_class, "Excluded");
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

    fn create_empty_clinvar() -> ClinVarAssessment {
        ClinVarAssessment {
            is_pathogenic: false,
            is_likely_pathogenic: false,
            is_benign: false,
            is_likely_benign: false,
            selected_entry: None,
            confidence_level: "none".to_string(),
            reason: "No ClinVar entries".to_string(),
        }
    }

    fn create_empty_predictive() -> PredictiveAssessment {
        PredictiveAssessment {
            suggests_pathogenic: false,
            contributing_scores: HashMap::new(),
            confidence: 0.0,
            support_count: 0,
            has_primate_ai_support: false,
        }
    }

    #[test]
    fn test_benign_filtering_with_pathogenic() {
        // Test that pathogenic variants are NOT filtered even with exclude_benign=true
        let variant = create_test_variant();
        let clinvar = ClinVarAssessment {
            is_pathogenic: true,
            is_likely_pathogenic: false,
            is_benign: true, // Has benign annotation too
            is_likely_benign: false,
            selected_entry: None,
            confidence_level: "high".to_string(),
            reason: "ClinVar pathogenic".to_string(),
        };
        let predictive = create_empty_predictive();

        let decision = make_filter_decision_with_config(&variant, &clinvar, &predictive, true);

        assert!(decision.should_include); // Should still be included
        assert_eq!(decision.pathogenicity_class, "Pathogenic");
        assert_eq!(decision.primary_evidence, "ClinVar");
    }

    #[test]
    fn test_benign_filtering_without_pathogenic() {
        // Test that benign-only variants ARE filtered when exclude_benign=true
        let variant = create_test_variant();
        let clinvar = ClinVarAssessment {
            is_pathogenic: false,
            is_likely_pathogenic: false,
            is_benign: true,
            is_likely_benign: false,
            selected_entry: None,
            confidence_level: "medium".to_string(),
            reason: "ClinVar benign".to_string(),
        };
        let predictive = create_empty_predictive();

        let decision = make_filter_decision_with_config(&variant, &clinvar, &predictive, true);

        assert!(!decision.should_include); // Should be excluded
        assert_eq!(decision.pathogenicity_class, "Excluded (Benign)");
        assert_eq!(decision.primary_evidence, "ClinVar");
    }

    #[test]
    fn test_benign_no_filtering_when_disabled() {
        // Test that benign variants are NOT filtered when exclude_benign=false
        let variant = create_test_variant();
        let clinvar = ClinVarAssessment {
            is_pathogenic: false,
            is_likely_pathogenic: false,
            is_benign: true,
            is_likely_benign: false,
            selected_entry: None,
            confidence_level: "medium".to_string(),
            reason: "ClinVar benign".to_string(),
        };
        let predictive = create_empty_predictive();

        let decision = make_filter_decision_with_config(&variant, &clinvar, &predictive, false);

        assert!(!decision.should_include); // Excluded for insufficient evidence, not benign
        assert_eq!(decision.pathogenicity_class, "Excluded");
        assert_eq!(decision.primary_evidence, "None");
    }
}
