use ahash::{AHashMap, HashMap};
use anyhow::Result;
use log::error;
use methylome::{Strand, find_motif_indices_in_contig, motif::Motif};
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

        let fwd_indices: Vec<usize> = find_motif_indices_in_contig(&contig_seq.as_str(), motif);
        let rev_indices: Vec<usize> =
            find_motif_indices_in_contig(&contig_seq.as_str(), &motif.reverse_complement());

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

        // // This is number of motif obervations with methylation data
        // let n_motif_obs = methylation_data.len() as u32;

        // let mean_read_cov = {
        //     let total_cov: u64 = methylation_data
        //         .iter()
        //         .map(|cov| cov.get_n_valid_cov() as u64)
        //         .sum();
        //     total_cov as f64 / methylation_data.len() as f64
        // };

        // let mut fractions: Vec<f64> = methylation_data
        //     .iter()
        //     .map(|cov| cov.fraction_modified())
        //     .collect();

        // fractions.sort_by(|a, b| a.partial_cmp(b).unwrap());
        // let median = if fractions.len() % 2 == 0 {
        //     let mid = fractions.len() / 2;
        //     (fractions[mid - 1] + fractions[mid]) / 2.0
        // } else {
        //     fractions[fractions.len() / 2]
        // };

        // local_results.push(MedianMotifMethylationDegree {
        //     contig: contig.id.clone(),
        //     motif: motif.clone(),
        //     median,
        //     mean_read_cov,
        //     n_motif_obs,
        //     motif_occurences_total,
        // })
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
    use csv::ReaderBuilder;
    use std::{
        fs::File,
        io::{BufReader, Write},
    };
    use tempfile::NamedTempFile;

    use crate::{
        models::genome_workspace::GenomeWorkspaceBuilder,
        services::domain::pileup_processor::parse_to_methylation_record,
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
            .add_contig(Contig::new(
                "contig_3".to_string(),
                "TGGACGATCCCGATC".to_string(),
            ))
            .unwrap();

        let file = File::open(pileup_file).unwrap();
        let reader = BufReader::new(file);
        let mut rdr = ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_reader(reader);

        for res in rdr.records() {
            let record = res.unwrap();

            let n_valid_cov_str = record.get(9).unwrap();
            let n_valid_cov = n_valid_cov_str.parse().unwrap();
            let meth_record =
                parse_to_methylation_record("contig_3".to_string(), &record, n_valid_cov, 0.8)
                    .unwrap();
            workspace_builder.add_record(meth_record.unwrap()).unwrap();
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
