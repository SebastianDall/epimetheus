use ahash::AHashMap;
use anyhow::anyhow;
use epimetheus_core::{
    models::{
        contig::Contig,
        genome_workspace::{GenomeWorkspace, GenomeWorkspaceBuilder},
        methylation::MethylationRecord,
        pileup::{PileupRecord, PileupRecordString},
    },
    services::traits::BatchLoader,
};
use log::{debug, warn};
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

pub struct SequentialBatchLoader<R: BufRead> {
    reader: R,
    assembly: AHashMap<String, Contig>,
    batch_size: usize,
    min_valid_read_coverage: u32,
    min_valid_cov_to_diff_fraction: f32,
    allow_mismatch: bool,

    current_contig_id: Option<String>,
    current_contig: Option<Contig>,
    pending_record: Option<Result<PileupRecordString, anyhow::Error>>,
    contigs_loaded_in_batch: usize,
}

impl<R: BufRead> SequentialBatchLoader<R> {
    pub fn new(
        reader: R,
        assembly: AHashMap<String, Contig>,
        batch_size: usize,
        min_valid_read_coverage: u32,
        min_valid_cov_to_diff_fraction: f32,
        allow_mismatch: bool,
    ) -> Self {
        let size = if batch_size == 0 {
            warn!("Batch size cannot be zero. Defaulting to 1.");
            1
        } else {
            batch_size
        };

        SequentialBatchLoader {
            reader,
            assembly,
            batch_size: size,
            min_valid_read_coverage,
            min_valid_cov_to_diff_fraction,
            allow_mismatch,
            current_contig_id: None,
            current_contig: None,
            pending_record: None,
            contigs_loaded_in_batch: 0,
        }
    }
}

impl BatchLoader<GenomeWorkspace> for SequentialBatchLoader<BufReader<File>> {
    fn next_batch(&mut self) -> Option<anyhow::Result<GenomeWorkspace>> {
        self.next()
    }

    fn new(
        reader: BufReader<File>,
        assembly: AHashMap<String, Contig>,
        batch_size: usize,
        min_valid_read_coverage: u32,
        min_valid_cov_to_diff_fraction: f32,
        allow_mismatch: bool,
    ) -> Self {
        Self::new(
            reader,
            assembly,
            batch_size,
            min_valid_read_coverage,
            min_valid_cov_to_diff_fraction,
            allow_mismatch,
        )
    }
}

impl<R: BufRead> Iterator for SequentialBatchLoader<R> {
    type Item = Result<GenomeWorkspace, anyhow::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut builder = GenomeWorkspaceBuilder::new();

        let record_iter = self
            .pending_record
            .take()
            .into_iter()
            .chain(std::iter::from_fn(|| {
                let mut line = String::new();
                match self.reader.read_line(&mut line) {
                    Ok(0) => None,
                    Ok(_) => {
                        let trimmed = line.trim_end().to_string();
                        if trimmed.is_empty() {
                            None
                        } else {
                            Some(Ok(PileupRecordString::new(trimmed)))
                        }
                    }
                    Err(e) => Some(Err(anyhow::Error::from(e))),
                }
            }));

        for record_result in record_iter {
            let record = match record_result {
                Ok(r) => r,
                Err(e) => return Some(Err(e)),
            };

            let record_for_pending = record.clone();
            let pileup_record = match PileupRecord::try_from(record) {
                Ok(p) => p,
                Err(e) => return Some(Err(e)),
            };

            let contig_id = pileup_record.contig.clone();

            if Some(&contig_id) != self.current_contig_id.as_ref() {
                debug!("Current contig id in line: {}", &contig_id);
                debug!(
                    "Current contig being added: {}",
                    self.current_contig
                        .as_ref()
                        .map(|c| c.id.to_string())
                        .unwrap_or("None".to_string())
                );

                match self.assembly.get(&contig_id) {
                    Some(found) => {
                        if let Some(old_contig) = self.current_contig.take() {
                            debug!("Adding contig to builder");
                            if let Err(e) = builder.add_contig(old_contig) {
                                return Some(Err(e));
                            }
                            self.contigs_loaded_in_batch += 1;
                            debug!(
                                "Contigs loaded in batch: {}. Batch size is: {}",
                                self.contigs_loaded_in_batch, self.batch_size
                            );

                            if self.contigs_loaded_in_batch == self.batch_size {
                                self.pending_record = Some(Ok(record_for_pending));
                                self.contigs_loaded_in_batch = 0;
                                return Some(Ok(builder.build()));
                            }
                        };

                        // Add the current contig to builder.
                        self.current_contig_id = Some(contig_id.clone());
                        self.current_contig = Some(found.clone());
                    }

                    // Return error if contig not found in assembly.
                    None if !self.allow_mismatch => {
                        return Some(Err(anyhow!("Contig '{}' not found in assembly", contig_id)));
                    }

                    // Skip records if mismatches are allowed
                    None => {
                        continue;
                    }
                }
            }
            let meth = match MethylationRecord::try_from_with_filters(
                pileup_record.clone(),
                self.min_valid_read_coverage,
                self.min_valid_cov_to_diff_fraction,
            ) {
                Ok(Some(m)) => m,
                Ok(None) => continue,
                Err(e) => return Some(Err(e)),
            };
            if let Some(ref mut c) = self.current_contig {
                if let Err(e) = c.add_methylation_record(meth) {
                    return Some(Err(e));
                }
            }
        }
        if let Some(last) = self.current_contig.take() {
            builder.add_contig(last).ok()?;
        }

        let workspace = builder.build();
        if workspace.is_empty() {
            None
        } else {
            Some(Ok(workspace))
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use epimetheus_core::models::methylation::MethylationCoverage;
    use std::{
        fs::File,
        io::{BufReader, Write},
    };
    use tempfile::NamedTempFile;

    #[test]
    fn test_batch_loading() -> anyhow::Result<()> {
        let mut pileup_file = NamedTempFile::new().unwrap();
        writeln!(
            pileup_file,
            "contig_3\t6\t1\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t8\t1\tm\t133\t+\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t12\t1\ta\t133\t+\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t7\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t13\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;

        let mut assembly = AHashMap::new();
        assembly.insert(
            "contig_3".to_string(),
            Contig::from_string("contig_3".to_string(), "TGGACGATCCCGATC".to_string()).unwrap(),
        );
        let file = File::open(pileup_file).unwrap();
        let reader = BufReader::new(file);

        let batch_loader = SequentialBatchLoader::new(reader, assembly, 1, 1, 0.8, false);

        for ws in batch_loader {
            let workspace = ws?.get_workspace();
            assert_eq!(workspace.len(), 1);

            assert_eq!(
                workspace
                    .get("contig_3")
                    .unwrap()
                    .get_methylated_positions(
                        &[6],
                        methylome::Strand::Positive,
                        methylome::ModType::SixMA
                    )
                    .into_iter()
                    .map(|(_pos, rec)| rec.unwrap().clone())
                    .collect::<Vec<MethylationCoverage>>(),
                vec![MethylationCoverage::new(15, 15, 0).unwrap()]
            );
        }

        Ok(())
    }

    #[test]
    fn test_multi_batches() -> anyhow::Result<()> {
        let mut pileup_file = NamedTempFile::new().unwrap();
        writeln!(
            pileup_file,
            "contig_3\t6\t1\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t8\t1\tm\t133\t+\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_4\t12\t1\ta\t133\t+\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_4\t7\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_4\t13\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;

        let mut assembly = AHashMap::new();
        assembly.insert(
            "contig_3".to_string(),
            Contig::from_string("contig_3".to_string(), "TGGACGATCCCGATC".to_string()).unwrap(),
        );
        assembly.insert(
            "contig_4".to_string(),
            Contig::from_string("contig_4".to_string(), "TGGACGATCCCGATC".to_string()).unwrap(),
        );
        let file = File::open(pileup_file).unwrap();
        let reader = BufReader::new(file);

        let batch_loader = SequentialBatchLoader::new(reader, assembly, 1, 1, 0.8, false);

        let mut num_batches = 0;
        for ws in batch_loader {
            num_batches += 1;

            if num_batches == 2 {
                let workspace = ws.unwrap().get_workspace();
                let contig_4 = workspace.get("contig_4").unwrap();
                assert_eq!(
                    contig_4
                        .get_methylated_positions(
                            &[12],
                            methylome::Strand::Positive,
                            methylome::ModType::SixMA
                        )
                        .into_iter()
                        .map(|(_pos, res)| res.unwrap().clone())
                        .collect::<Vec<MethylationCoverage>>(),
                    vec![MethylationCoverage::new(5, 20, 0).unwrap()]
                );
            }
        }

        assert_eq!(num_batches, 2);

        Ok(())
    }

    #[test]
    fn test_contig_missing_error() -> anyhow::Result<()> {
        let mut pileup_file = NamedTempFile::new().unwrap();
        writeln!(
            pileup_file,
            "contig_3\t6\t1\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t8\t1\tm\t133\t+\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_4\t12\t1\ta\t133\t+\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_4\t7\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_4\t13\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;

        let assembly = AHashMap::new();
        let file = File::open(pileup_file).unwrap();
        let reader = BufReader::new(file);

        let batch_loader = SequentialBatchLoader::new(reader, assembly, 2, 1, 0.8, false);

        for ws in batch_loader {
            assert!(ws.is_err());
        }

        Ok(())
    }

    #[test]
    fn test_contig_mismatch_exits() -> anyhow::Result<()> {
        let mut pileup_file = NamedTempFile::new().unwrap();
        writeln!(
            pileup_file,
            "contig_3\t6\t1\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t8\t1\tm\t133\t+\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_5\t12\t1\ta\t133\t+\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?; // This contig does not exist in the assembly
        writeln!(
            pileup_file,
            "contig_4\t7\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_4\t13\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;

        let mut assembly = AHashMap::new();
        assembly.insert(
            "contig_3".to_string(),
            Contig::from_string("contig_3".to_string(), "TGGACGATCCCGATC".to_string()).unwrap(),
        );
        assembly.insert(
            "contig_4".to_string(),
            Contig::from_string("contig_4".to_string(), "TGGACGATCCCGATC".to_string()).unwrap(),
        );
        let file = File::open(pileup_file).unwrap();
        let reader = BufReader::new(file);

        let batch_loader = SequentialBatchLoader::new(reader, assembly, 2, 1, 0.8, false);

        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            for ws in batch_loader {
                let _workspace = ws.unwrap(); // This should fail on the third record
            }
        }));

        assert!(
            result.is_err(),
            "Expected the program to panic due to missing contig, but it did not."
        );

        Ok(())
    }

    #[test]
    fn test_contig_mismatch_exits_but_is_allowed() -> anyhow::Result<()> {
        let mut pileup_file = NamedTempFile::new().unwrap();
        writeln!(
            pileup_file,
            "contig_3\t6\t1\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_3\t8\t1\tm\t133\t+\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_5\t12\t1\ta\t133\t+\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?; // This contig does not exist in the assembly
        writeln!(
            pileup_file,
            "contig_4\t7\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t20\t123\t0\t0\t6\t0\t0"
        )?;
        writeln!(
            pileup_file,
            "contig_4\t13\t1\ta\t133\t-\t0\t1\t255,0,0\t20\t0.00\t5\t123\t0\t0\t6\t0\t0"
        )?;

        let mut assembly = AHashMap::new();
        assembly.insert(
            "contig_3".to_string(),
            Contig::from_string("contig_3".to_string(), "TGGACGATCCCGATC".to_string()).unwrap(),
        );
        assembly.insert(
            "contig_4".to_string(),
            Contig::from_string("contig_4".to_string(), "TGGACGATCCCGATC".to_string()).unwrap(),
        );
        let file = File::open(pileup_file).unwrap();
        let reader = BufReader::new(file);

        let batch_loader = SequentialBatchLoader::new(reader, assembly, 3, 1, 0.8, true);

        for ws in batch_loader {
            assert_eq!(ws.unwrap().get_workspace().len(), 2);
        }

        Ok(())
    }
}
