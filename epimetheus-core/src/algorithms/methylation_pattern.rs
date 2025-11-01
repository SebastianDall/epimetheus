use ahash::{AHashMap, HashMap};
use anyhow::Result;
use log::error;
use methylome::{Strand, find_motif_indices_in_sequence, motif::Motif};
use rayon::prelude::*;

use crate::models::{
    contig::{Contig, ContigId, Position as ContigPosition},
    genome_workspace::GenomeWorkspace,
    methylation::{MethylationCoverage, MotifMethylationPositions},
};

pub fn calculate_contig_read_methylation_single(
    contig: &Contig,
    motifs: Vec<Motif>,
) -> Result<MotifMethylationPositions> {
    let contig_seq = &contig.sequence;

    let mut all_methylation_data = AHashMap::new();

    for motif in motifs.iter() {
        let mod_type = motif.mod_type;

        let fwd_indices: Vec<usize> = find_motif_indices_in_sequence(&contig_seq, motif);
        let rev_indices: Vec<usize> =
            find_motif_indices_in_sequence(&contig_seq, &motif.reverse_complement());

        if fwd_indices.is_empty() && rev_indices.is_empty() {
            continue;
        }

        // This is the actual number of motifs in the contig
        // let motif_occurences_total = fwd_indices.len() as u32 + rev_indices.len() as u32;

        let fwd_methylation =
            contig.get_methylated_positions(&fwd_indices, methylome::Strand::Positive, mod_type);
        let rev_methylation =
            contig.get_methylated_positions(&rev_indices, methylome::Strand::Negative, mod_type);

        let methylation_data_fwd: HashMap<
            (ContigId, Motif, ContigPosition, Strand),
            MethylationCoverage,
        > = fwd_methylation
            .into_iter()
            .filter_map(|(pos, maybe_cov)| {
                maybe_cov.map(|meth| {
                    (
                        (contig.id.clone(), motif.clone(), pos, Strand::Positive),
                        meth.clone(),
                    )
                })
            })
            .collect();

        let methylation_data_rev: HashMap<
            (ContigId, Motif, ContigPosition, Strand),
            MethylationCoverage,
        > = rev_methylation
            .into_iter()
            .filter_map(|(pos, maybe_cov)| {
                maybe_cov.map(|meth| {
                    (
                        (contig.id.clone(), motif.clone(), pos, Strand::Negative),
                        meth.clone(),
                    )
                })
            })
            .collect();

        if methylation_data_rev.is_empty() & methylation_data_fwd.is_empty() {
            continue;
        }

        all_methylation_data.extend(methylation_data_fwd);
        all_methylation_data.extend(methylation_data_rev);
    }

    Ok(MotifMethylationPositions {
        methylation: all_methylation_data,
    })
}

pub fn calculate_contig_read_methylation_pattern(
    contigs: GenomeWorkspace,
    motifs: Vec<Motif>,
    num_threads: usize,
) -> Result<MotifMethylationPositions> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Could not initialize threadpool");

    let mut combined_contig_motif_methylation = AHashMap::new();
    let results: Vec<MotifMethylationPositions> = contigs
        .get_workspace()
        .par_iter()
        .map(|(contig_id, contig)| {
            calculate_contig_read_methylation_single(contig, motifs.clone()).unwrap_or_else(|e| {
                error!("Error processing contig {}: {}", contig_id, e);
                MotifMethylationPositions::new(AHashMap::new())
            })
        })
        .collect();

    for res in results {
        combined_contig_motif_methylation.extend(res.methylation);
    }

    Ok(MotifMethylationPositions::new(
        combined_contig_motif_methylation,
    ))
}

#[cfg(test)]
mod tests {
    use std::{
        fs::File,
        io::{BufRead, BufReader, Write},
    };
    use tempfile::NamedTempFile;

    use crate::models::{
        genome_workspace::GenomeWorkspaceBuilder,
        methylation::MethylationRecord,
        pileup::{PileupRecord, PileupRecordString},
    };

    use super::*;

    #[test]
    fn test_calculate_methylation() -> Result<()> {
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

        let mut workspace_builder = GenomeWorkspaceBuilder::new();

        // Add a mock contig to the workspace
        workspace_builder
            .add_contig(
                Contig::from_string("contig_3".to_string(), "TGGACGATCCCGATC".to_string()).unwrap(),
            )
            .unwrap();

        let file = File::open(pileup_file).unwrap();
        let reader = BufReader::new(file);

        for res in reader.lines() {
            let record = res.unwrap();
            let pileup_record = PileupRecord::try_from(PileupRecordString::new(record)).unwrap();
            let meth_record = MethylationRecord::try_from_with_filters(pileup_record, 1, 0.8)?;
            if let Some(meth) = meth_record {
                workspace_builder.add_record(meth).unwrap();
            }
        }

        let workspace = workspace_builder.build();

        let motifs = vec![
            Motif::new("GATC", "a", 1).unwrap(),
            Motif::new("GATC", "m", 3).unwrap(),
            Motif::new("GATC", "21839", 3).unwrap(),
        ];
        let contig_methylation_pattern =
            calculate_contig_read_methylation_pattern(workspace, motifs, 1).unwrap();

        let expected_median_result = vec![0.625, 1.0];
        let mut meth_result_median: Vec<f64> = contig_methylation_pattern
            .to_median_degrees()
            .iter()
            .map(|res| res.median)
            .collect();
        meth_result_median.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(meth_result_median, expected_median_result);

        let expected_weighted_mean_result = vec![0.6, 1.0];
        let mut meth_result_weighted_mean: Vec<f64> = contig_methylation_pattern
            .to_weighted_mean_degress()
            .iter()
            .map(|res| res.w_mean)
            .collect();
        meth_result_weighted_mean.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(meth_result_weighted_mean, expected_weighted_mean_result);

        let expected_mean_read_cov = vec![18.75, 20.0];
        let mut meth_result: Vec<f64> = contig_methylation_pattern
            .to_median_degrees()
            .iter()
            .map(|res| res.mean_read_cov)
            .collect();
        meth_result.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(meth_result, expected_mean_read_cov);

        let expected_n_motif_obs = vec![1, 4];
        let mut meth_result: Vec<u32> = contig_methylation_pattern
            .to_median_degrees()
            .iter()
            .map(|res| res.n_motif_obs)
            .collect();
        meth_result.sort_by(|a, b| a.cmp(b));
        assert_eq!(meth_result, expected_n_motif_obs);

        Ok(())
    }
}
