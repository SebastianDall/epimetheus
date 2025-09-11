use std::{path::Path, sync::Mutex};

use ahash::AHashMap;
use anyhow::Result;
use csv::StringRecord;
use epimetheus_support::bgzip::reader::PileupReader;
use rayon::prelude::*;

use crate::{
    data::{GenomeWorkspace, GenomeWorkspaceBuilder, contig::Contig},
    extract_methylation_pattern::{parse_to_methylation_record, reader::BatchReader},
};

pub struct ParallelBatchLoader {
    readers: Vec<PileupReader>,
    assembly: AHashMap<String, Contig>,
    batch_size: usize,
    min_valid_read_coverage: u32,
    min_valid_cov_to_diff_fraction: f32,
    allow_mismatch: bool,

    // Iterator fields
    processed_contigs: Option<Vec<String>>,
}

impl ParallelBatchLoader {
    pub fn new(
        file: &Path,
        assembly: AHashMap<String, Contig>,
        batch_size: usize,
        min_valid_read_coverage: u32,
        min_valid_cov_to_diff_fraction: f32,
        allow_mismatch: bool,
        threads: usize,
    ) -> Result<Self> {
        let readers: Result<Vec<_>> = (0..threads)
            .map(|_| PileupReader::from_path(file))
            .collect();

        Ok(Self {
            readers: readers?,
            assembly,
            batch_size,
            min_valid_read_coverage,
            min_valid_cov_to_diff_fraction,
            allow_mismatch,
            processed_contigs: None,
        })
    }
}

impl Iterator for ParallelBatchLoader {
    type Item = Result<GenomeWorkspace, anyhow::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut builder = GenomeWorkspaceBuilder::new();

        let contigs_in_pileup = self.readers[0].available_contigs();

        let contigs_to_be_processed = match &self.processed_contigs {
            Some(processed_contigs) => contigs_in_pileup
                .iter()
                .filter(|c| !processed_contigs.contains(c))
                .cloned()
                .collect::<Vec<String>>(),
            None => contigs_in_pileup,
        };

        let batch: Vec<String> = contigs_to_be_processed
            .into_iter()
            .filter(|contig_id| {
                if self.allow_mismatch {
                    self.assembly.contains_key(contig_id)
                } else {
                    true
                }
            })
            .take(self.batch_size)
            .collect();

        if batch.is_empty() {
            return None;
        }

        let readers_mutex: Vec<Mutex<&mut PileupReader>> =
            self.readers.iter_mut().map(|r| Mutex::new(r)).collect();

        let batch_results: Result<Vec<Contig>, anyhow::Error> = batch
            .into_par_iter()
            .enumerate()
            .map(|(i, contig_id)| {
                let reader_index = i % readers_mutex.len();
                let mut reader = readers_mutex[reader_index].lock().unwrap();

                let assembly_contig = self.assembly.get(&contig_id).expect("Contig should exist in assembly after filtering. Consider using allow assembly pileup mismatch.");

                process_contig(
                    &mut **reader,
                    assembly_contig,
                    self.min_valid_read_coverage,
                    self.min_valid_cov_to_diff_fraction,
                )
            })
            .collect();

        let mut processed_contigs = Vec::new();
        match batch_results {
            Ok(res) => {
                for contig in res {
                    let contig_id = contig.id.clone();
                    processed_contigs.push(contig_id.clone());
                    builder.add_contig(contig).expect(&format!(
                        "Error adding contig '{}' to builder. This should be infallible..",
                        contig_id
                    ));
                }
            }
            Err(e) => return Some(Err(e)),
        }

        if let Some(ref mut list) = self.processed_contigs {
            list.extend(processed_contigs);
        } else {
            self.processed_contigs = Some(processed_contigs);
        }

        let workspace = builder.build();
        if workspace.is_empty() {
            None
        } else {
            Some(Ok(workspace))
        }
    }
}

impl BatchReader for ParallelBatchLoader {
    fn next_batch(&mut self) -> Option<Result<GenomeWorkspace>> {
        self.next()
    }
}

pub fn process_contig(
    reader: &mut PileupReader,
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
