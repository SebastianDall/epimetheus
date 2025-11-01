use anyhow::Result;
use methylome::{Motif, find_motif_indices_in_sequence, read::Read};
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;

pub fn extract_read_methylation_pattern(
    reads: Vec<Read>,
    motifs: Vec<Motif>,
    min_meth_qual: u8,
    threads: usize,
) -> Result<DataFrame> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Could not initialize threadpool");

    // Process reads in batches to control memory usage
    const BATCH_SIZE: usize = 1000;
    let batches: Vec<_> = reads.chunks(BATCH_SIZE).collect();

    // Process batches in parallel
    let results: Vec<HashMap<(String, String, String, u32), (u32, u32)>> = batches
        .into_par_iter()
        .map(|batch| {
            let mut batch_data = HashMap::new();

            for read in batch {
                let sequence = read.get_sequence();
                let modifications = read.get_modifications();

                for motif in &motifs {
                    // Find all motif positions in this read
                    let indices = find_motif_indices_in_sequence(sequence, motif);

                    if !indices.is_empty() {
                        let motif_sequence = motif
                            .sequence
                            .iter()
                            .map(|b| b.to_string())
                            .collect::<String>();

                        let key = (
                            read.get_name().clone(),
                            motif_sequence,
                            motif.mod_type.to_pileup_code().to_string(),
                            motif.mod_position as u32,
                        );

                        let mut n_modified = 0;
                        let n_motif_obs = indices.len() as u32;

                        // Count modifications that meet quality threshold
                        for &position in &indices {
                            if let Some(meth_base) = modifications.0.get(&position) {
                                if meth_base.quality.0 >= min_meth_qual {
                                    n_modified += 1;
                                }
                            }
                        }

                        // Aggregate counts within this batch
                        let entry = batch_data.entry(key).or_insert((0, 0));
                        entry.0 += n_modified;
                        entry.1 += n_motif_obs;
                    }
                }
            }
            batch_data
        })
        .collect();

    // Merge results from all batches
    let mut aggregated_data: HashMap<(String, String, String, u32), (u32, u32)> = HashMap::new();
    for batch_result in results {
        for (key, (n_modified, n_motif_obs)) in batch_result {
            let entry = aggregated_data.entry(key).or_insert((0, 0));
            entry.0 += n_modified;
            entry.1 += n_motif_obs;
        }
    }

    // Convert aggregated data to vectors for DataFrame
    let mut read_ids = Vec::with_capacity(aggregated_data.len());
    let mut motif_sequences = Vec::with_capacity(aggregated_data.len());
    let mut mod_types = Vec::with_capacity(aggregated_data.len());
    let mut mod_positions = Vec::with_capacity(aggregated_data.len());
    let mut n_modified_vec = Vec::with_capacity(aggregated_data.len());
    let mut n_motif_obs_vec = Vec::with_capacity(aggregated_data.len());

    for ((read_id, motif_seq, mod_type, mod_pos), (n_modified, n_motif_obs)) in aggregated_data {
        read_ids.push(read_id);
        motif_sequences.push(motif_seq);
        mod_types.push(mod_type);
        mod_positions.push(mod_pos);
        n_modified_vec.push(n_modified);
        n_motif_obs_vec.push(n_motif_obs);
    }

    // Create DataFrame
    let df = df! [
        "read_id" => read_ids,
        "motif" => motif_sequences,
        "mod_type" => mod_types,
        "mod_position" => mod_positions,
        "n_modified" => n_modified_vec,
        "n_motif_obs" => n_motif_obs_vec,
    ]?;

    Ok(df)
}
