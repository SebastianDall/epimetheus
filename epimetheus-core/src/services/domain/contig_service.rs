use anyhow::Result;

use crate::models::{contig::Contig, methylation::MethylationRecord};

pub fn populate_contig_with_methylation(
    contig: &Contig,
    records: Vec<MethylationRecord>,
) -> Result<Contig> {
    let mut contig = contig.clone();

    for rec in records {
        contig.add_methylation_record(rec)?;
    }
    Ok(contig)
}
