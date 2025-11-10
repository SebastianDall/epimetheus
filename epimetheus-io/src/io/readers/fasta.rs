use ahash::AHashMap;
use anyhow::{Context, anyhow};
use epimetheus_core::models::contig::Contig;
use methylome::sequence::Sequence;
use seq_io::fasta::{Reader as FxReader, Record};
use std::path::Path;

use crate::io::traits::FastaReader;

pub struct Reader;

impl FastaReader for Reader {
    fn read_fasta(
        path: &Path,
        contig_filter: Option<Vec<String>>,
    ) -> anyhow::Result<AHashMap<String, Contig>> {
        let mut fasta_reader = FxReader::from_path(&path)
            .with_context(|| format!("Failed to open FASTA at: {:?}", path))?;

        let mut contigs = AHashMap::new();

        while let Some(record_result) = fasta_reader.next() {
            let record = record_result.with_context(|| "Error reading record from FASTA file.")?;

            let id = record
                .id()
                .map(String::from)
                .with_context(|| "Error extracting record ID")?;

            if let Some(ref contig_filter) = contig_filter {
                if !contig_filter.contains(&id) {
                    continue;
                }
            }

            let seq = Sequence::from_u8(record.seq())
                .map_err(|e| anyhow!("Could not parse contig '{}': {}", id, e.to_string()))?;

            contigs.insert(id.clone(), Contig::new(id, seq));
        }
        Ok(contigs)
    }
}
