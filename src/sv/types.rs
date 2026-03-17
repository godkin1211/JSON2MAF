/// Core data types for SV (Structural Variant) processing.
///
/// SVs differ fundamentally from SNVs:
/// - Positions are ranges (start..svEnd), not points
/// - Read support comes from split reads and paired-end reads
/// - Predictive scores (REVEL, DANN, PrimateAI) do not apply
/// - HGVSc/HGVSp only available for sequence-resolved SVs, not symbolic (<DEL>, <INS>)
use crate::TranscriptAnnotation;

// ============================================================================
// SV Classification
// ============================================================================

#[derive(Debug, Clone, PartialEq)]
pub enum SVType {
    Del,
    Ins,
}

impl SVType {
    pub fn as_str(&self) -> &'static str {
        match self {
            SVType::Del => "DEL",
            SVType::Ins => "INS",
        }
    }
}

// ============================================================================
// ClinGen Annotation
// ============================================================================

/// ClinGen dosage sensitivity annotation, present at the position level for SVs.
/// Multiple entries may exist; we select the most clinically relevant.
#[derive(Debug, Clone)]
pub struct ClinGenEntry {
    pub id: String,
    pub variant_type: String,
    pub clinical_interpretation: String,
    pub phenotypes: Vec<String>,
}

// ============================================================================
// Read Support
// ============================================================================

/// Read-level evidence for an SV call.
///
/// `split_read_counts[0]` = reads spanning the breakpoint supporting the alt allele.
/// `paired_end_read_counts[0]` = discordant read pairs supporting the alt allele.
/// Index [1] in each is the reference-supporting count.
#[derive(Debug, Clone, Default)]
pub struct SVReadSupport {
    pub split_read_alt: i64,
    pub split_read_ref: i64,
    pub paired_end_alt: i64,
    pub paired_end_ref: i64,
}

impl SVReadSupport {
    pub fn total_alt(&self) -> i64 {
        self.split_read_alt + self.paired_end_alt
    }

    pub fn total_ref(&self) -> i64 {
        self.split_read_ref + self.paired_end_ref
    }

    /// VAF computed from all alt-supporting reads over total depth estimate.
    pub fn vaf(&self) -> Option<f64> {
        let total = self.total_alt() + self.total_ref();
        if total > 0 {
            Some(self.total_alt() as f64 / total as f64)
        } else {
            None
        }
    }
}

// ============================================================================
// Parsed SV Position
// ============================================================================

/// A fully parsed SV position, ready for output conversion.
#[derive(Debug, Clone)]
pub struct SVPosition {
    pub chromosome: String,
    pub start: i64,
    pub end_pos: i64,
    pub sv_type: SVType,
    /// Length in base pairs of the SV event.
    pub sv_length: i64,
    /// True for symbolic alleles (<DEL>, <INS>); false for sequence-resolved.
    /// Symbolic SVs generally lack HGVSc/HGVSp.
    pub is_symbolic: bool,
    pub reference_allele: String,
    pub alternate_allele: String,
    pub filters: Vec<String>,
    pub read_support: SVReadSupport,
    /// Transcript annotations from variants[0]. May be empty for some symbolic SVs.
    pub transcripts: Vec<TranscriptAnnotation>,
    pub clingen: Vec<ClinGenEntry>,
    pub sample_name: String,
}

// ============================================================================
// Output Record (TSV row)
// ============================================================================

/// One row in the SV output TSV file.
#[derive(Debug, Clone)]
pub struct SVRecord {
    pub hugo_symbol: String,
    pub chromosome: String,
    pub start_position: i64,
    pub end_position: i64,
    pub sv_type: String,
    pub sv_length: i64,
    /// Most severe consequence across all affected transcripts.
    pub variant_classification: String,
    pub hgvsc: String,
    pub hgvsp: String,
    pub transcript_id: String,
    pub split_read_alt: i64,
    pub split_read_ref: i64,
    pub paired_end_alt: i64,
    pub paired_end_ref: i64,
    pub total_alt_support: i64,
    pub total_ref_support: i64,
    pub vaf: String,
    pub filters: String,
    pub tumor_sample_barcode: String,
    pub clingen_id: String,
    pub clingen_interpretation: String,
    pub clingen_phenotypes: String,
}
