use anyhow::Result;
use methylome::{Motif, find_motif_indices_in_sequence, read::Read};
use polars::prelude::*;

pub fn extract_read_methylation_pattern(
    reads: Vec<Read>,
    motifs: Vec<Motif>,
    // _min_meth_quality: MethQual,
) -> Result<DataFrame> {
    let mut read_ids = Vec::new();
    let mut motif_sequences = Vec::new();
    let mut mod_positions = Vec::new();
    let mut mod_types = Vec::new();
    let mut read_positions = Vec::new();
    let mut meth_quals = Vec::new();

    for read in reads {
        let sequence = read.get_sequence();
        let modifications = read.get_modifications();

        for motif in &motifs {
            // Find all motif positions in this read
            let indices = find_motif_indices_in_sequence(sequence, motif);

            for &position in &indices {
                read_ids.push(read.get_name().clone());
                motif_sequences.push(
                    motif
                        .sequence
                        .iter()
                        .map(|b| b.to_string())
                        .collect::<String>(),
                );
                mod_positions.push(motif.mod_position as u32);
                mod_types.push(motif.mod_type.to_pileup_code());
                read_positions.push(position as u32);

                // Check if this position is methylated
                if let Some(meth_base) = modifications.0.get(&position) {
                    meth_quals.push(Some(meth_base.quality.0 as u32));
                } else {
                    meth_quals.push(None);
                }
            }
        }
    }

    // Create DataFrame
    let df = df! [
        "read_id" => read_ids,
        "motif" => motif_sequences,
        "mod_type" => mod_types,
        "mod_position" => mod_positions,
        "read_position" => read_positions,
        "meth_qual" => meth_quals,
    ]?;

    Ok(df)
}
