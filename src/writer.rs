use crate::types::MAFRecord;
use anyhow::{Context, Result};
use csv::Writer;
use std::fs::File;
use std::path::Path;

pub struct MAFWriter {
    writer: Writer<File>,
    records_written: usize,
}

impl MAFWriter {
    pub fn new(output_path: &str) -> Result<Self> {
        let file = File::create(output_path)
            .with_context(|| format!("Failed to create output file: {}", output_path))?;

        // Create writer with headers enabled
        // The csv crate will automatically write headers from the serde field names
        // on the first call to serialize()
        let writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_writer(file);

        Ok(Self {
            writer,
            records_written: 0,
        })
    }

    pub fn write_record(&mut self, record: &MAFRecord) -> Result<()> {
        self.writer
            .serialize(record)
            .context("Failed to write MAF record")?;
        self.records_written += 1;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush().context("Failed to flush writer")?;
        Ok(())
    }

    pub fn records_written(&self) -> usize {
        self.records_written
    }
}

pub fn merge_maf_files(input_files: &[String], output_path: &str) -> Result<usize> {
    let mut output = MAFWriter::new(output_path)?;
    let mut total_records = 0;

    for input_file in input_files {
        if !Path::new(input_file).exists() {
            continue;
        }

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(input_file)
            .with_context(|| format!("Failed to open input file: {}", input_file))?;

        for result in reader.deserialize() {
            let record: MAFRecord = result.context("Failed to deserialize MAF record")?;
            output.write_record(&record)?;
            total_records += 1;
        }
    }

    output.flush()?;
    Ok(total_records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_write_maf() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let output_path = temp_dir.path().join("test.maf");
        let output_str = output_path.to_str().unwrap();

        let mut writer = MAFWriter::new(output_str)?;

        let record = MAFRecord {
            hugo_symbol: "BRAF".to_string(),
            chromosome: "chr7".to_string(),
            start_position: 140453136,
            end_position: 140453136,
            strand: "+".to_string(),
            variant_classification: "Missense_Mutation".to_string(),
            variant_type: "SNP".to_string(),
            reference_allele: "A".to_string(),
            tumor_seq_allele1: "A".to_string(),
            tumor_seq_allele2: "T".to_string(),
            tumor_sample_barcode: "SAMPLE1".to_string(),
            hgvsc: "c.1799T>A".to_string(),
            hgvsp: "p.Val600Glu".to_string(),
            hgvsp_short: "p.V600E".to_string(),
            transcript_id: "NM_004333.4".to_string(),
            exon: "15/18".to_string(),
            consequence: "missense_variant".to_string(),
            impact: "MODERATE".to_string(),
            codons: "Gtg/Gag".to_string(),
            amino_acids: "V/E".to_string(),
            cdna_position: "1799/2301".to_string(),
            cds_position: "1799/2301".to_string(),
            protein_position: "600/766".to_string(),
            dbsnp_rs: "rs113488022".to_string(),
            dbsnp_val_status: "".to_string(),
            cosmic_id: "COSM476".to_string(),
            clinvar_id: "RCV000123456".to_string(),
            clinvar_review_status: "reviewed by expert panel".to_string(),
            clinvar_significance: "Pathogenic".to_string(),
            clinvar_disease: "Cancer".to_string(),
            primate_ai_score: "0.85".to_string(),
            dann_score: "0.99".to_string(),
            revel_score: "0.92".to_string(),
            gnomad_af: "0.0001".to_string(),
            gnomad_eas_af: "0.0".to_string(),
            depth: "100".to_string(),
            vaf: "0.45".to_string(),
        };

        writer.write_record(&record)?;
        writer.flush()?;

        assert_eq!(writer.records_written(), 1);
        assert!(output_path.exists());

        Ok(())
    }
}
