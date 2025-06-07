use anyhow::{anyhow, Context, Result};
use rayon::prelude::*;
use log::info;
use methylome::{IupacBase, Motif};
use std::{
    collections::HashSet,
    io::{BufWriter, Write},
    path::Path,
};

pub mod args;

pub use args::MotifClusteringArgs;

use crate::{processing::create_motifs, utils::create_output_file};

fn pick_victim(m1: &Motif, m2: &Motif) -> Motif {
    let len1 = m1.sequence_to_string().len();
    let len2 = m2.sequence_to_string().len();
    if len1 > len2 {
        m1.clone()
    } else if len1 < len2 {
        m2.clone()
    } else if m1.possible_dna_sequences().len() > m2.possible_dna_sequences().len() {
        m2.clone()
    } else {
        m1.clone()
    }
}

fn collapse_child_motifs(motifs: &[Motif]) -> Vec<Motif> {
    let n = motifs.len();

    // 1) in parallel, scan all (i,j) pairs and collect your “victims”
    let victims: Vec<Motif> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            // for each i, scan j = i+1..n in parallel
            (i + 1..n)
                .into_par_iter()
                .filter_map(move |j| {
                    let m1 = &motifs[i];
                    let m2 = &motifs[j];
                    if m1.is_child_motif(m2) || m2.is_child_motif(m1) {
                        // pick the shorter/less‐possible one
                        let victim = pick_victim(&m1, &m2);
                        Some(victim)
                    } else {
                        None
                    }
                })
        })
        .collect();

    // 2) turn your victims into a HashSet for O(1) lookups
    let remove_set: HashSet<Motif> = victims.into_iter().collect();

    // 3) in parallel, keep only those not in remove_set
    motifs
        .par_iter()
        .filter(|m| !remove_set.contains(*m))
        .cloned()
        .collect()
}


#[allow(dead_code)]
fn collapse_motifs(motifs: &Vec<Motif>) -> Result<Motif> {
    let first_motif = motifs[0].clone();
    let n_bases = first_motif.sequence.len();

    for m in motifs {
        if m.sequence.len() != n_bases {
            return Err(anyhow!("Not all motifs have the same length"));
        } else if m.mod_type != first_motif.mod_type {
            return Err(anyhow!("Not all motifs have the same modification"));
        } else if m.mod_position != first_motif.mod_position {
            return Err(anyhow!(
                "Motifs does not have the same mod_position. Cannot create final motif: {:#?}",
                motifs
            ));
        }
    }

    let mut sequence = Vec::with_capacity(n_bases);
    for i in 0..n_bases {
        let mut nucs = HashSet::new();
        for motif in motifs {
            for possible_nuc in motif.sequence[i].to_possible_nucleotides() {
                nucs.insert(possible_nuc);
            }
        }
        let unified_base = IupacBase::from_nucleotides(&nucs)?;
        sequence.push(unified_base);
    }

    let seq = sequence
        .iter()
        .map(IupacBase::to_string)
        .collect::<Vec<_>>()
        .join("");

    let final_motif = Motif::new(
        seq.as_str(),
        first_motif.mod_type.to_pileup_code(),
        first_motif.mod_position,
    )?;

    Ok(final_motif)
}

pub fn motif_clustering(args: MotifClusteringArgs) -> Result<()> {
    let outpath = Path::new(&args.output);

    create_output_file(outpath)?;

    let motifs = match args.motifs {
        Some(motifs) => {
            info!("Motifs loaded");
            create_motifs(motifs).context("Failed to parse motifs")?
        }
        _ => {
            anyhow::bail!("No motifs found");
        }
    };

    let motifs_with_no_childs = collapse_child_motifs(&motifs);

    // let mut uf = cluster_motifs(&motifs_with_no_childs, true);
    // let motif_clusters = group_motifs_by_set(&mut uf, &motifs_with_no_childs);

    // // Within cluster find best candidate motif
    // // Should be the smallest
    // // In case several have the same lenght
    // // they should be collapsed
    // let mut motif_cluster_representatives = HashMap::new();

    // for (_cluster, motifs_in_cluster) in motif_clusters {
    //     let min_motif = motifs_in_cluster
    //         .iter()
    //         .map(|m| m.sequence.len())
    //         .min()
    //         .unwrap();
    //     let smallest_motifs = motifs_in_cluster
    //         .iter()
    //         .cloned()
    //         .filter(|m| m.sequence.len() == min_motif)
    //         .collect::<Vec<_>>();

    //     if smallest_motifs.len() > 1 {
    //         let mut rep_cluster = cluster_motifs(&smallest_motifs, true);
    //         let rep_motif_clusters = group_motifs_by_set(&mut rep_cluster, &smallest_motifs);

    //         for (_rep_cluster, rep_motifs_in_cluster) in rep_motif_clusters {
    //             let rep_motif = collapse_motifs(&rep_motifs_in_cluster)?;
    //             motif_cluster_representatives.insert(rep_motif, motifs_in_cluster.clone());
    //         }
    //     } else {
    //         let rep_motif = smallest_motifs[0].clone();
    //         motif_cluster_representatives.insert(rep_motif, motifs_in_cluster);
    //     }
    // }

    let outfile = std::fs::File::create(outpath)
        .with_context(|| format!("Failed to create file at: {:?}", outpath))?;
    let mut writer = BufWriter::new(outfile);

    // writeln!(
    //     writer,
    //     "motif_representative\tmod_type_representative\tmod_position_representative\tmotif\tmod_type\tmod_position"
    // )?;

    // for (rep, motifs) in motif_cluster_representatives {
    //     let rep_motif_sequence = rep.sequence_to_string();
    //     let rep_mod_type_str = rep.mod_type.to_pileup_code();
    //     let rep_mod_position = rep.mod_position;

    //     for motif in motifs {
    //         let motif_sequence = motif.sequence_to_string();
    //         let mod_type_str = motif.mod_type.to_pileup_code();
    //         let mod_position = motif.mod_position;

    //         writeln!(
    //             writer,
    //             "{}\t{}\t{}\t{}\t{}\t{}",
    //             rep_motif_sequence,
    //             rep_mod_type_str,
    //             rep_mod_position,
    //             motif_sequence,
    //             mod_type_str,
    //             mod_position,
    //         )?;
    //         writer.flush()?;
    //     }
    // }


    for m in motifs_with_no_childs {
        writeln!(writer, "{}\t{}\t{}", m.sequence_to_string(), m.mod_type.to_pileup_code(), m.mod_position)?;
    }

    
    Ok(())
}

#[cfg(test)]
mod tests {
    use methylome::Motif;

    use super::*;

    #[test]
    fn test_collapse_child_motifs() {
        let m1 = Motif::new("GATC", "m", 3).unwrap();
        let m2 = Motif::new("GGATC", "m", 4).unwrap();
        let m3 = Motif::new("GTTCT", "m", 3).unwrap();
        let m4 = Motif::new("GATCC", "m", 3).unwrap();
        let m5 = Motif::new("GATC", "a", 1).unwrap();

        let motifs = vec![m1.clone(), m2.clone(), m3.clone(), m4.clone(), m5.clone()];

        let motifs_to_keep = collapse_child_motifs(&motifs);

        assert_eq!(motifs_to_keep.len(), 3);
        assert_eq!(motifs_to_keep[0], m1.clone());
        assert_eq!(motifs_to_keep[1], m3.clone());
        assert_eq!(motifs_to_keep[2], m5.clone());
    }

    

}
