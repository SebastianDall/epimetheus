use std::path::Path;

use ahash::AHashMap;
use anyhow::Result;
use epimetheus_core::{
    models::{contig::Contig, genome_workspace::GenomeWorkspace, pileup::PileupRecord},
    services::traits::BatchLoader,
};

use crate::io::traits::PileupReader;

pub fn load_pileup_records_for_contig<R: PileupReader>(
    pileup_path: &Path,
    contig_id: &str,
) -> anyhow::Result<Vec<PileupRecord>> {
    let mut reader = R::from_path(pileup_path)?;
    let pileup_record_strings = reader.query_contig(contig_id)?;
    pileup_record_strings
        .into_iter()
        .map(PileupRecord::try_from)
        .collect::<anyhow::Result<Vec<PileupRecord>>>()
}

pub fn process_batches_from_loader<L: BatchLoader<GenomeWorkspace>>(
    loader: &mut L,
) -> impl Iterator<Item = Result<AHashMap<String, Contig>>> + '_ {
    std::iter::from_fn(move || match loader.next_batch() {
        Some(Ok(workspace)) => {
            let contigs = workspace.get_workspace();
            Some(Ok(contigs))
        }
        Some(Err(e)) => Some(Err(e)),
        None => None,
    })
}
