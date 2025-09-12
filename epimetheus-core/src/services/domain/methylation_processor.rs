use anyhow::Result;
use csv::StringRecord;

use crate::{
    models::contig::Contig,
    services::{domain::pileup_processor::parse_to_methylation_record, traits::PileupReader},
};

pub fn process_contig<R: PileupReader>(
    reader: &mut R,
    assembly_contig: &Contig,
    min_valid_read_coverage: u32,
    min_valid_cov_to_diff_fraction: f32,
) -> Result<Contig> {
    let records = reader.query_contig(&assembly_contig.id)?;
    let mut contig = assembly_contig.clone();

    for record in records {
        let pileup_line = StringRecord::from(record.0.split('\t').collect::<Vec<&str>>());
        let methylation_record = parse_to_methylation_record(
            contig.id.clone(),
            &pileup_line,
            min_valid_read_coverage,
            min_valid_cov_to_diff_fraction,
        )?;

        if let Some(meth) = methylation_record {
            contig.add_methylation_record(meth)?;
        } else {
            continue;
        }
    }

    Ok(contig)
}
