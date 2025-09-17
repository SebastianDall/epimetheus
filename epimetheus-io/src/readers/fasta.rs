use ahash::AHashMap;
use anyhow::Context;
use epimetheus_core::{models::contig::Contig, services::traits::FastaReader};
use seq_io::fasta::{Reader as FxReader, Record};
use std::path::Path;

pub struct Reader;

impl FastaReader for Reader {
    fn read_fasta(path: &Path) -> anyhow::Result<AHashMap<String, Contig>> {
        let mut fasta_reader = FxReader::from_path(&path)
            .with_context(|| format!("Failed to open FASTA at: {:?}", path))?;

        let mut contigs = AHashMap::new();

        while let Some(record_result) = fasta_reader.next() {
            let record = record_result.with_context(|| "Error reading record from FASTA file.")?;

            let id = record
                .id()
                .map(String::from)
                .with_context(|| "Error extracting record ID")?;

            let seq = String::from_utf8(record.owned_seq())
                .with_context(|| format!("Invalid UTF8 character in FASTA record: '{}'", id))?
                .to_string();

            contigs.insert(id.clone(), Contig::new(id, seq));
        }
        Ok(contigs)
    }
}
