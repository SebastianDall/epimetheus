use anyhow::Result;
use epimetheus_core::models::pileup::PileupRecord;

use crate::io::traits::PileupReader;

pub fn query_pileup<R: PileupReader>(
    reader: &mut R,
    contigs: &[String],
) -> Result<Vec<PileupRecord>> {
    let mut all_records = Vec::new();

    for c in contigs {
        let records = reader.query_contig(&c)?;

        for rec in records {
            let pileup_rec = PileupRecord::try_from(rec)?;
            all_records.push(pileup_rec);
        }
    }
    Ok(all_records)
}
