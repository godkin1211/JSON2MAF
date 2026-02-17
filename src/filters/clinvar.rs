use crate::types::*;

pub fn assess_clinvar_pathogenicity(entries: &[ClinVarEntry]) -> ClinVarAssessment {
    // If no ClinVar annotations, return negative result
    if entries.is_empty() {
        return ClinVarAssessment {
            is_pathogenic: false,
            is_likely_pathogenic: false,
            is_benign: false,
            is_likely_benign: false,
            selected_entry: None,
            confidence_level: "none".to_string(),
            reason: "No ClinVar entries available".to_string(),
        };
    }

    // Filter Pathogenic and Likely pathogenic entries
    let pathogenic_entries: Vec<_> = entries
        .iter()
        .filter(|e| is_pathogenic_entry(e))
        .collect();

    // Filter Benign and Likely benign entries
    let benign_entries: Vec<_> = entries
        .iter()
        .filter(|e| is_benign_entry(e))
        .collect();

    // Check for benign classification first
    let (is_benign, is_likely_benign) = if !benign_entries.is_empty() {
        let selected_benign = resolve_conflicting_entries(&benign_entries);
        let sig_lower = selected_benign
            .clinical_significance
            .join(", ")
            .to_lowercase();
        let parts: Vec<&str> = sig_lower.split(&['/', ',', ';'][..]).collect();

        let has_standalone_benign = parts
            .iter()
            .any(|p| p.trim().contains("benign") && !p.trim().contains("likely"));
        let has_likely_benign = parts
            .iter()
            .any(|p| p.trim().contains("likely") && p.trim().contains("benign"));

        (has_standalone_benign, !has_standalone_benign && has_likely_benign)
    } else {
        (false, false)
    };

    // If no pathogenic entries
    if pathogenic_entries.is_empty() {
        return ClinVarAssessment {
            is_pathogenic: false,
            is_likely_pathogenic: false,
            is_benign,
            is_likely_benign,
            selected_entry: None,
            confidence_level: "none".to_string(),
            reason: "No pathogenic or likely pathogenic ClinVar entries".to_string(),
        };
    }

    // Select best entry
    let selected = resolve_conflicting_entries(&pathogenic_entries).clone();

    // Determine pathogenicity level
    let sig_lower = selected
        .clinical_significance
        .join(", ")
        .to_lowercase();

    // Check if contains standalone "Pathogenic" (not just "Likely pathogenic")
    let parts: Vec<&str> = sig_lower.split(&['/', ',', ';'][..]).collect();
    let has_standalone_pathogenic = parts
        .iter()
        .any(|p| p.trim().contains("pathogenic") && !p.trim().contains("likely"));
    let has_likely_pathogenic = parts
        .iter()
        .any(|p| p.trim().contains("likely") && p.trim().contains("pathogenic"));

    let is_path = has_standalone_pathogenic;
    let is_likely_path = !has_standalone_pathogenic && has_likely_pathogenic;

    // Determine confidence level
    let confidence = get_confidence_level(selected.review_status.as_deref().unwrap_or(""));

    // Build assessment reason
    let reason = build_assessment_reason(&selected, pathogenic_entries.len());

    ClinVarAssessment {
        is_pathogenic: is_path,
        is_likely_pathogenic: is_likely_path,
        is_benign,
        is_likely_benign,
        selected_entry: Some(selected),
        confidence_level: confidence,
        reason,
    }
}

fn is_pathogenic_entry(entry: &ClinVarEntry) -> bool {
    let sig_lower = entry.clinical_significance.join(", ").to_lowercase();
    sig_lower.contains("pathogenic")
        && !sig_lower.contains("benign")
        && !sig_lower.contains("uncertain")
}

fn is_benign_entry(entry: &ClinVarEntry) -> bool {
    let sig_lower = entry.clinical_significance.join(", ").to_lowercase();
    sig_lower.contains("benign")
        && !sig_lower.contains("pathogenic")
        && !sig_lower.contains("uncertain")
}

pub fn get_review_status_priority(status: &str) -> i32 {
    let status_lower = status.to_lowercase();

    if status_lower.contains("practice guideline") {
        1
    } else if status_lower.contains("reviewed by expert panel") {
        2
    } else if status_lower.contains("multiple submitters") && status_lower.contains("no conflict")
    {
        3
    } else if status_lower.contains("conflicting") {
        4
    } else if status_lower.contains("single submitter") {
        5
    } else if status_lower.contains("no assertion criteria provided") {
        6
    } else if status_lower.contains("no assertion provided") {
        7
    } else {
        8 // Unknown status, lowest priority
    }
}

pub fn is_cancer_related(diseases: &[String]) -> bool {
    let cancer_keywords = [
        "cancer",
        "carcinoma",
        "tumor",
        "tumour",
        "malignant",
        "neoplasm",
        "lymphoma",
        "leukemia",
        "leukaemia",
        "sarcoma",
        "melanoma",
        "glioma",
        "blastoma",
        "myeloma",
        "adenocarcinoma",
    ];

    diseases.iter().any(|disease| {
        let disease_lower = disease.to_lowercase();
        cancer_keywords
            .iter()
            .any(|keyword| disease_lower.contains(keyword))
    })
}

pub fn resolve_conflicting_entries<'a>(entries: &'a [&'a ClinVarEntry]) -> &'a ClinVarEntry {
    if entries.len() == 1 {
        return entries[0];
    }

    // Sort by priority
    let mut sorted = entries.to_vec();
    sorted.sort_by(|a, b| {
        // First compare by review status priority (lower is better)
        let pa = get_review_status_priority(a.review_status.as_deref().unwrap_or(""));
        let pb = get_review_status_priority(b.review_status.as_deref().unwrap_or(""));
        if pa != pb {
            return pa.cmp(&pb);
        }

        // Then compare by cancer-related (cancer-related comes first)
        let ca = is_cancer_related(&a.phenotypes);
        let cb = is_cancer_related(&b.phenotypes);
        if ca != cb {
            return cb.cmp(&ca); // Reverse order for bool
        }

        // Finally compare by last_evaluated (newer comes first)
        let da = a.last_evaluated.as_deref().unwrap_or("");
        let db = b.last_evaluated.as_deref().unwrap_or("");
        db.cmp(da)
    });

    sorted[0]
}

fn get_confidence_level(review_status: &str) -> String {
    let priority = get_review_status_priority(review_status);

    if priority <= 2 {
        "high".to_string()
    } else if priority <= 4 {
        "medium".to_string()
    } else {
        "low".to_string()
    }
}

fn build_assessment_reason(entry: &ClinVarEntry, total_entries: usize) -> String {
    let mut parts = Vec::new();

    // ClinVar classification
    if !entry.clinical_significance.is_empty() {
        parts.push(format!("ClinVar: {}", entry.clinical_significance.join(", ")));
    }

    // Review status
    if let Some(status) = &entry.review_status {
        parts.push(format!("Review: {}", status));
    }

    // If multiple entries
    if total_entries > 1 {
        parts.push(format!("Selected from {} entries", total_entries));
    }

    // Cancer-related
    if is_cancer_related(&entry.phenotypes) {
        parts.push("Cancer-related disease".to_string());
    }

    parts.join("; ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_pathogenic_entry() {
        let entry = ClinVarEntry {
            id: Some("RCV001".to_string()),
            allele_id: None,
            clinical_significance: vec!["Pathogenic".to_string()],
            review_status: None,
            phenotypes: vec![],
            last_evaluated: None,
        };
        assert!(is_pathogenic_entry(&entry));
    }

    #[test]
    fn test_is_cancer_related() {
        let diseases = vec!["Breast cancer".to_string()];
        assert!(is_cancer_related(&diseases));

        let diseases = vec!["Diabetes".to_string()];
        assert!(!is_cancer_related(&diseases));
    }

    #[test]
    fn test_review_status_priority() {
        assert_eq!(get_review_status_priority("practice guideline"), 1);
        assert_eq!(get_review_status_priority("reviewed by expert panel"), 2);
        assert_eq!(
            get_review_status_priority("criteria provided, multiple submitters, no conflicts"),
            3
        );
    }
}
