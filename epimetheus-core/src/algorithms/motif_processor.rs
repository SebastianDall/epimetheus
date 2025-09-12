use anyhow::{Result, anyhow};
use methylome::{IupacBase, Motif};
use rayon::prelude::*;
use std::collections::HashSet;

#[allow(dead_code)]
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

pub fn collapse_child_motifs(motifs: &[Motif]) -> Vec<Motif> {
    let n = motifs.len();

    // 1) in parallel, scan all (i,j) pairs and collect your “victims”
    let victims: Vec<Motif> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            // for each i, scan j = i+1..n in parallel
            (i + 1..n).into_par_iter().filter_map(move |j| {
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

#[cfg(test)]
mod tests {
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
