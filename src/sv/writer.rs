/// Tab-separated output writer for SV records.
use anyhow::Result;
use std::fs::File;
use std::io::{BufWriter, Write};

use super::types::SVRecord;

pub const SV_TSV_HEADERS: &[&str] = &[
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "SV_Type",
    "SV_Length",
    "Variant_Classification",
    "HGVSc",
    "HGVSp",
    "Transcript_ID",
    "Split_Read_Alt",
    "Split_Read_Ref",
    "Paired_End_Alt",
    "Paired_End_Ref",
    "Total_Alt_Support",
    "Total_Ref_Support",
    "VAF",
    "Filters",
    "Tumor_Sample_Barcode",
    "ClinGen_ID",
    "ClinGen_Interpretation",
    "ClinGen_Phenotypes",
];

pub struct SVWriter {
    inner: BufWriter<File>,
}

impl SVWriter {
    pub fn new(path: &str) -> Result<Self> {
        let file = File::create(path)?;
        let mut inner = BufWriter::new(file);
        writeln!(inner, "{}", SV_TSV_HEADERS.join("\t"))?;
        Ok(Self { inner })
    }

    pub fn write_record(&mut self, r: &SVRecord) -> Result<()> {
        writeln!(
            self.inner,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.hugo_symbol,
            r.chromosome,
            r.start_position,
            r.end_position,
            r.sv_type,
            r.sv_length,
            r.variant_classification,
            r.hgvsc,
            r.hgvsp,
            r.transcript_id,
            r.split_read_alt,
            r.split_read_ref,
            r.paired_end_alt,
            r.paired_end_ref,
            r.total_alt_support,
            r.total_ref_support,
            r.vaf,
            r.filters,
            r.tumor_sample_barcode,
            r.clingen_id,
            r.clingen_interpretation,
            r.clingen_phenotypes,
        )?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.inner.flush().map_err(Into::into)
    }
}
