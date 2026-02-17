use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// ============================================================================
// Filter Configuration
// ============================================================================

#[derive(Debug, Clone)]
pub struct FilterConfig {
    // Quality filtering parameters
    pub min_total_depth: i32,
    pub min_variant_frequency: f64,

    // Population frequency filtering parameters
    pub max_eas_af: f64,

    // Predictive score thresholds
    pub min_revel_score: f64,
    pub min_primate_ai_score: f64,
    pub min_dann_score: f64,

    // ClinVar filtering options
    pub exclude_benign: bool,
}

impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            min_total_depth: 30,
            min_variant_frequency: 0.03,
            max_eas_af: 0.01,
            min_revel_score: 0.75,
            min_primate_ai_score: 0.8,
            min_dann_score: 0.96,
            exclude_benign: false,
        }
    }
}

impl FilterConfig {
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.min_total_depth < 1 {
            anyhow::bail!("min_total_depth must be at least 1, got {}", self.min_total_depth);
        }

        if !(0.0..=1.0).contains(&self.min_variant_frequency) {
            anyhow::bail!("min_variant_frequency must be between 0 and 1, got {}", self.min_variant_frequency);
        }

        if !(0.0..=1.0).contains(&self.max_eas_af) {
            anyhow::bail!("max_eas_af must be between 0 and 1, got {}", self.max_eas_af);
        }

        if !(0.0..=1.0).contains(&self.min_revel_score) {
            anyhow::bail!("min_revel_score must be between 0 and 1, got {}", self.min_revel_score);
        }

        if !(0.0..=1.0).contains(&self.min_primate_ai_score) {
            anyhow::bail!("min_primate_ai_score must be between 0 and 1, got {}", self.min_primate_ai_score);
        }

        if !(0.0..=1.0).contains(&self.min_dann_score) {
            anyhow::bail!("min_dann_score must be between 0 and 1, got {}", self.min_dann_score);
        }

        Ok(())
    }
}

// ============================================================================
// Nirvana JSON Data Structures
// ============================================================================

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct DataSource {
    pub name: String,
    pub version: String,
    pub description: Option<String>,
    pub release_date: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct NirvanaHeader {
    pub annotator: String,
    pub creation_time: String,
    pub genome_assembly: String,
    pub schema_version: i32,
    pub data_sources: Vec<DataSource>,
    pub samples: Vec<String>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ClinVarEntry {
    pub id: Option<String>,
    #[serde(rename = "variantId")]
    pub allele_id: Option<String>,
    #[serde(rename = "significance", default)]
    pub clinical_significance: Vec<String>,
    #[serde(rename = "reviewStatus")]
    pub review_status: Option<String>,
    #[serde(default)]
    pub phenotypes: Vec<String>,
    pub last_evaluated: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct TranscriptAnnotation {
    #[serde(rename = "transcript")]
    pub id: Option<String>,
    pub source: Option<String>,
    pub hgnc: Option<String>,
    #[serde(default)]
    pub consequence: Vec<String>,
    pub impact: Option<String>,
    pub amino_acids: Option<String>,
    pub cdna_pos: Option<String>,
    pub cds_pos: Option<String>,
    pub exons: Option<String>,
    pub codons: Option<String>,
    pub protein_pos: Option<String>,
    pub hgvsc: Option<String>,
    pub hgvsp: Option<String>,
    #[serde(rename = "isCanonical")]
    pub is_canonical: Option<bool>,
    #[serde(rename = "isManeSelect")]
    pub is_mane_select: Option<bool>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PopulationFrequency {
    #[serde(rename = "population")]
    pub source: String,
    #[serde(rename = "allAf")]
    pub all_af: Option<f64>,
    #[serde(rename = "easAf")]
    pub eas_af: Option<f64>,
    #[serde(rename = "afrAf")]
    pub afr_af: Option<f64>,
    #[serde(rename = "amrAf")]
    pub amr_af: Option<f64>,
    #[serde(rename = "eurAf")]
    pub eur_af: Option<f64>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct CosmicEntry {
    pub id: Option<String>,
    pub gene: Option<String>,
    #[serde(rename = "mutationType")]
    pub mutation_type: Option<String>,
    pub count: Option<i32>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Sample {
    pub total_depth: Option<i32>,
    pub variant_frequencies: Option<Vec<f64>>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PrimateAIEntry {
    pub hgnc: Option<String>,
    pub score_percentile: Option<f64>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PrimateAI3DEntry {
    pub score: Option<f64>,
    pub score_percentile: Option<f64>,
    pub classification: Option<String>,
    pub ensembl_transcript_id: Option<String>,
    pub ref_seq_transcript_id: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct RevelScore {
    pub score: Option<f64>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Variant {
    pub variant_type: String,
    #[serde(default)]
    pub transcripts: Vec<TranscriptAnnotation>,
    #[serde(default)]
    pub clinvar: Vec<ClinVarEntry>,
    #[serde(default)]
    pub cosmic: Vec<CosmicEntry>,
    #[serde(default)]
    pub dbsnp: Vec<String>,
    #[serde(rename = "primateAI-3D", default)]
    pub primate_ai_3d: Vec<PrimateAI3DEntry>,
    #[serde(rename = "primateAI", default)]
    pub primate_ai: Vec<PrimateAIEntry>,
    #[serde(rename = "dannScore")]
    pub dann_score: Option<f64>,
    #[serde(rename = "revel")]
    pub revel_score: Option<RevelScore>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Position {
    pub chromosome: String,
    pub position: i32,
    #[serde(rename = "refAllele")]
    pub reference_allele: String,
    #[serde(rename = "altAlleles")]
    pub alternate_alleles: Vec<String>,
    #[serde(default)]
    pub filters: Vec<String>,
    #[serde(default)]
    pub samples: Vec<Sample>,
    #[serde(default)]
    pub variants: Vec<Variant>,
}

#[derive(Debug, Clone)]
pub struct VariantPosition {
    pub chromosome: String,
    pub start: i32,
    pub end_pos: i32,
    pub reference_allele: String,
    pub alternate_allele: String,
    pub variant_type: String,

    // Filter status
    pub filters: Vec<String>,

    // Sample information
    pub total_depth: Option<i32>,
    pub variant_frequencies: Option<Vec<f64>>,

    // Annotation information
    pub transcripts: Vec<TranscriptAnnotation>,
    pub clinvar: Vec<ClinVarEntry>,
    pub cosmic: Vec<CosmicEntry>,
    pub population_frequencies: Vec<PopulationFrequency>,

    // Predictive scores
    pub primate_ai_3d: Option<f64>,
    pub primate_ai: Option<f64>,
    pub dann_score: Option<f64>,
    pub revel_score: Option<f64>,

    // dbSNP
    pub dbsnp_ids: Vec<String>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct NirvanaData {
    pub header: NirvanaHeader,
    pub positions: Vec<Position>,
}

// ============================================================================
// Filter Assessment Results
// ============================================================================

#[derive(Debug, Clone)]
pub struct QualityFilterResult {
    pub passes_quality: bool,
    pub failure_reason: Option<String>,
    pub depth: Option<i32>,
    pub variant_frequency: Option<f64>,
    pub eas_allele_frequency: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct ClinVarAssessment {
    pub is_pathogenic: bool,
    pub is_likely_pathogenic: bool,
    pub is_benign: bool,
    pub is_likely_benign: bool,
    pub selected_entry: Option<ClinVarEntry>,
    pub confidence_level: String,
    pub reason: String,
}

#[derive(Debug, Clone)]
pub struct PredictiveAssessment {
    pub suggests_pathogenic: bool,
    pub contributing_scores: HashMap<String, f64>,
    pub confidence: f64,
    pub support_count: usize,
    pub has_primate_ai_support: bool,
}

#[derive(Debug, Clone)]
pub struct FilterDecision {
    pub should_include: bool,
    pub pathogenicity_class: String,
    pub primary_evidence: String,
    pub justification: String,
}

// ============================================================================
// MAF Format
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MAFRecord {
    #[serde(rename = "Hugo_Symbol")]
    pub hugo_symbol: String,
    #[serde(rename = "Chromosome")]
    pub chromosome: String,
    #[serde(rename = "Start_Position")]
    pub start_position: i32,
    #[serde(rename = "End_Position")]
    pub end_position: i32,
    #[serde(rename = "Strand")]
    pub strand: String,
    #[serde(rename = "Variant_Classification")]
    pub variant_classification: String,
    #[serde(rename = "Variant_Type")]
    pub variant_type: String,
    #[serde(rename = "Reference_Allele")]
    pub reference_allele: String,
    #[serde(rename = "Tumor_Seq_Allele1")]
    pub tumor_seq_allele1: String,
    #[serde(rename = "Tumor_Seq_Allele2")]
    pub tumor_seq_allele2: String,
    #[serde(rename = "Tumor_Sample_Barcode")]
    pub tumor_sample_barcode: String,
    #[serde(rename = "HGVSc")]
    pub hgvsc: String,
    #[serde(rename = "HGVSp")]
    pub hgvsp: String,
    #[serde(rename = "HGVSp_Short")]
    pub hgvsp_short: String,
    #[serde(rename = "Transcript_ID")]
    pub transcript_id: String,
    #[serde(rename = "Exon")]
    pub exon: String,
    #[serde(rename = "Consequence")]
    pub consequence: String,
    #[serde(rename = "IMPACT")]
    pub impact: String,
    #[serde(rename = "Codons")]
    pub codons: String,
    #[serde(rename = "Amino_Acids")]
    pub amino_acids: String,
    #[serde(rename = "cDNA_position")]
    pub cdna_position: String,
    #[serde(rename = "CDS_position")]
    pub cds_position: String,
    #[serde(rename = "Protein_position")]
    pub protein_position: String,
    #[serde(rename = "dbSNP_RS")]
    pub dbsnp_rs: String,
    #[serde(rename = "dbSNP_Val_Status")]
    pub dbsnp_val_status: String,
    #[serde(rename = "COSMIC_ID")]
    pub cosmic_id: String,
    #[serde(rename = "ClinVar_ID")]
    pub clinvar_id: String,
    #[serde(rename = "ClinVar_Review_Status")]
    pub clinvar_review_status: String,
    #[serde(rename = "ClinVar_Significance")]
    pub clinvar_significance: String,
    #[serde(rename = "ClinVar_Disease")]
    pub clinvar_disease: String,
    #[serde(rename = "PrimateAI_Score")]
    pub primate_ai_score: String,
    #[serde(rename = "DANN_Score")]
    pub dann_score: String,
    #[serde(rename = "REVEL_Score")]
    pub revel_score: String,
    #[serde(rename = "gnomAD_AF")]
    pub gnomad_af: String,
    #[serde(rename = "gnomAD_EAS_AF")]
    pub gnomad_eas_af: String,
    #[serde(rename = "Depth")]
    pub depth: String,
    #[serde(rename = "VAF")]
    pub vaf: String,
}

// ============================================================================
// Statistics
// ============================================================================

#[derive(Debug, Clone, Default)]
pub struct FilterStats {
    pub passed_quality: usize,
    pub failed_depth: usize,
    pub failed_vaf: usize,
    pub failed_af: usize,
    pub clinvar_pathogenic: usize,
    pub clinvar_likely: usize,
    pub predictive_likely: usize,
    pub primate_ai_only: usize,
    pub multi_score: usize,
    pub excluded_benign: usize,
    pub included: usize,
    pub excluded: usize,
}

impl FilterStats {
    pub fn merge(&mut self, other: &FilterStats) {
        self.passed_quality += other.passed_quality;
        self.failed_depth += other.failed_depth;
        self.failed_vaf += other.failed_vaf;
        self.failed_af += other.failed_af;
        self.clinvar_pathogenic += other.clinvar_pathogenic;
        self.clinvar_likely += other.clinvar_likely;
        self.predictive_likely += other.predictive_likely;
        self.primate_ai_only += other.primate_ai_only;
        self.multi_score += other.multi_score;
        self.excluded_benign += other.excluded_benign;
        self.included += other.included;
        self.excluded += other.excluded;
    }
}
