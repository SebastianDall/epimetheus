use ahash::AHashMap;
use anyhow::anyhow;
use log::{debug, warn};
use methylome::{ModType, Strand};
use rayon::prelude::*;
use std::io::{BufRead, BufReader};

use crate::data::{
    contig::Contig, methylation::MethylationCoverage, GenomeWorkspace, GenomeWorkspaceBuilder,
    MethylationRecord,
};

pub struct BatchLoader<R> {
    reader: BufReader<R>,
    assembly: AHashMap<String, Contig>,
    batch_size: usize,
    min_valid_read_coverage: u32,

    current_contig_id: Option<String>,
    current_contig: Option<Contig>,
    contigs_loaded_in_batch: usize,
    done: bool,
}

impl<R: std::io::Read> BatchLoader<R> {
    pub fn new(
        inner: R,
        assembly: AHashMap<String, Contig>,
        batch_size: usize,
        min_valid_read_coverage: u32,
    ) -> Self {
        let size = if batch_size == 0 {
            warn!("Batch size cannot be zero. Defaulting to 1.");
            1
        } else {
            batch_size
        };

        BatchLoader {
            reader: BufReader::new(inner),
            assembly,
            batch_size: size,
            min_valid_read_coverage,
            current_contig_id: None,
            current_contig: None,
            contigs_loaded_in_batch: 0,
            done: false,
        }
    }

    fn read_pileup_chunk(&mut self, chunk_size: usize) -> std::io::Result<Option<Vec<String>>> {
        let mut lines = Vec::with_capacity(chunk_size);
        for _ in 0..chunk_size {
            let mut buf = String::new();
            let n = self.reader.read_line(&mut buf)?;
            if n == 0 {
                break;
            }
            lines.push(buf);
        }

        if lines.is_empty() {
            Ok(None)
        } else {
            Ok(Some(lines))
        }
    }
}

impl<R: BufRead> Iterator for BatchLoader<R> {
    type Item = Result<GenomeWorkspace, anyhow::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let mut builder = GenomeWorkspaceBuilder::new();

        loop {
            if self.contigs_loaded_in_batch >= self.batch_size {
                self.contigs_loaded_in_batch = 0;
                return Some(Ok(builder.build()));
            }

            let maybe_pileup_chunk = match self.read_pileup_chunk(100_000) {
                Ok(mc) => mc,
                Err(e) => return Some(Err(e.into())),
            };

            let pileup_chunk = match maybe_pileup_chunk {
                Some(c) => c,
                None => {
                    if let Some(last_contig) = self.current_contig.take() {
                        if let Err(e) = builder.add_contig(last_contig) {
                            return Some(Err(e));
                        }
                    }

                    let workspace = builder.build();

                    if workspace.is_empty() {
                        self.done = true;
                        return None;
                    } else {
                        self.done = true;
                        return Some(Ok(workspace));
                    }
                }
            };

            let parsed_records: Result<Vec<_>, anyhow::Error> = pileup_chunk
                .par_iter()
                .map(|line| {
                    let cols: Vec<&str> = line.trim_end().split('\t').collect();
                    if cols.len() != 18 {
                        return Err(anyhow!("Not enough columns in pileup chunk"));
                    }

                    let contig_id = cols[0].to_string();
                    let n_valid_cov: u32 = cols[9]
                        .parse()
                        .map_err(|_| anyhow!("Invalid n_valid_coverage number"))?;
                    let position: usize = cols[1]
                        .parse()
                        .map_err(|_| anyhow!("Invalid mod_position field"))?;
                    let mod_type: ModType = cols[3].parse()?;
                    let strand: Strand = cols[5].parse()?;
                    let n_modified: u32 = cols[11]
                        .parse()
                        .map_err(|_| anyhow!("Invalid n_modified field"))?;

                    let methylation = MethylationCoverage::new(n_modified, n_valid_cov)?;

                    let methylation_record =
                        MethylationRecord::new(contig_id, position, strand, mod_type, methylation);

                    Ok(methylation_record)
                })
                .filter_map(|meth_rec| match meth_rec {
                    Ok(r)
                        if r.get_methylation_coverage().get_n_valid_cov()
                            >= self.min_valid_read_coverage =>
                    {
                        Some(Ok(r))
                    }
                    Ok(_) => None,
                    Err(e) => Some(Err(e)),
                })
                .collect::<Result<Vec<_>, anyhow::Error>>();

            let records = match parsed_records {
                Ok(r) => r,
                Err(e) => return Some(Err(e)),
            };

            for meth_rec in records {
                let contig_id = meth_rec.get_contig_id();
                // If the contig id in MethylationRecord does not match current contig_id being build -> switch contig
                if Some(&contig_id) != self.current_contig_id.as_ref() {
                    debug!("Current contig id in line: {}", &contig_id);
                    debug!(
                        "Current contig being added: {}",
                        self.current_contig
                            .as_ref()
                            .map(|c| c.id.to_string())
                            .unwrap_or("None".to_string())
                    );

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
                    };
                    self.current_contig_id = Some(contig_id.clone());

                    let contig_key = self.current_contig_id.as_ref().unwrap();
                    let new_contig = match self.assembly.get(contig_key) {
                        Some(c) => c.clone(),
                        None => {
                            return Some(Err(anyhow!(
                                "Contig '{}' not found in assembly",
                                contig_key
                            )))
                        }
                    };
                    self.current_contig = Some(new_contig);
                }
                if let Some(ref mut c) = self.current_contig {
                    if let Err(e) = c.add_methylation_record(meth_rec) {
                        return Some(Err(e));
                    }
                }
            }
        }
    }
}

// for record_result in self.reader.records() {
//     let record = match record_result.context("Failed to read pileup record") {
//         Ok(r) => r,
//         Err(e) => return Some(Err(e)),
//     };

//     let n_valid_cov: u32 = match record.get(9) {
//         Some(f) => match f.parse() {
//             Ok(v) => v,
//             Err(_) => return Some(Err(anyhow!("Invalid coverage number"))),
//         },
//         None => return Some(Err(anyhow!("Invalid coverage field"))),
//     };

//     if n_valid_cov < self.min_valid_read_coverage {
//         continue;
//     }

//     let contig_id = match record.get(0) {
//         Some(id) => id.to_string(),
//         None => return Some(Err(anyhow!("Missing contig id field"))),
//     };

//     if Some(&contig_id) != self.current_contig_id.as_ref() {
//         debug!("Current contig id in line: {}", &contig_id);
//         debug!(
//             "Current contig being added: {}",
//             self.current_contig
//                 .as_ref()
//                 .map(|c| c.id.to_string())
//                 .unwrap_or("None".to_string())
//         );

//         // Add the current contig to builder.
//         if let Some(old_contig) = self.current_contig.take() {
//             debug!("Adding contig to builder");
//             if let Err(e) = builder.add_contig(old_contig) {
//                 return Some(Err(e));
//             }
//             self.contigs_loaded_in_batch += 1;
//             debug!(
//                 "Contigs loaded in batch: {}. Batch size is: {}",
//                 self.contigs_loaded_in_batch, self.batch_size
//             );
//         };

//         self.current_contig_id = Some(contig_id.clone());

//         let contig_key = self.current_contig_id.as_ref().unwrap();
//         let new_contig = match self.assembly.get(contig_key) {
//             Some(c) => c.clone(),
//             None => {
//                 return Some(Err(anyhow!(
//                     "Contig '{}' not found in assembly",
//                     contig_key
//                 )))
//             }
//         };
//         self.current_contig = Some(new_contig);
//     }
//     let meth = match parse_to_methylation_record(contig_id, n_valid_cov, &record) {
//         Ok(m) => m,
//         Err(e) => return Some(Err(e)),
//     };
//     if let Some(ref mut c) = self.current_contig {
//         if let Err(e) = c.add_methylation_record(meth) {
//             return Some(Err(e));
//         }
//     }

//     if self.contigs_loaded_in_batch == self.batch_size {
//         self.contigs_loaded_in_batch = 0;
//         return Some(Ok(builder.build()));
//     }
// }
// if let Some(last) = self.current_contig.take() {
//     builder.add_contig(last).ok()?;
// }

// let workspace = builder.build();
// if workspace.is_empty() {
//     None
// } else {
//     Some(Ok(workspace))
// }
#[cfg(test)]
mod tests {
    use crate::data::methylation::MethylationCoverage;

    use super::*;
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
            Contig::new("contig_3".to_string(), "TGGACGATCCCGATC".to_string()),
        );
        let file = File::open(pileup_file).unwrap();
        let reader = BufReader::new(file);

        let batch_loader = BatchLoader::new(reader, assembly, 1, 1);

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
                    .map(|rec| rec.unwrap().clone())
                    .collect::<Vec<MethylationCoverage>>(),
                vec![MethylationCoverage::new(15, 15).unwrap()]
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
            Contig::new("contig_3".to_string(), "TGGACGATCCCGATC".to_string()),
        );
        assembly.insert(
            "contig_4".to_string(),
            Contig::new("contig_4".to_string(), "TGGACGATCCCGATC".to_string()),
        );
        let file = File::open(pileup_file).unwrap();
        let reader = BufReader::new(file);

        let batch_loader = BatchLoader::new(reader, assembly, 1, 1);

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
                        .map(|res| res.unwrap().clone())
                        .collect::<Vec<MethylationCoverage>>(),
                    vec![MethylationCoverage::new(5, 20).unwrap()]
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

        let batch_loader = BatchLoader::new(reader, assembly, 2, 1);

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
            Contig::new("contig_3".to_string(), "TGGACGATCCCGATC".to_string()),
        );
        assembly.insert(
            "contig_4".to_string(),
            Contig::new("contig_4".to_string(), "TGGACGATCCCGATC".to_string()),
        );
        let file = File::open(pileup_file).unwrap();
        let reader = BufReader::new(file);

        let batch_loader = BatchLoader::new(reader, assembly, 2, 1);

        let result = std::panic::catch_unwind(|| {
            for ws in batch_loader {
                let _workspace = ws.unwrap(); // This should fail on the third record
            }
        });

        assert!(
            result.is_err(),
            "Expected the program to panic due to missing contig, but it did not."
        );

        Ok(())
    }
}
