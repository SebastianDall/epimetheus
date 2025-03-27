use anyhow::{anyhow, Context, Result};
use log::info;
use methylome::{IupacBase, Motif};
use std::{
    collections::{HashMap, HashSet},
    io::{BufWriter, Write},
    path::Path,
};

pub mod args;

pub use args::MotifClusteringArgs;

use crate::{processing::create_motifs, utils::create_output_file};

struct UnionFind {
    pub parent: Vec<usize>,
    pub rank: Vec<usize>,
}

impl UnionFind {
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    pub fn find(&mut self, i: usize) -> usize {
        if self.parent[i] != i {
            self.parent[i] = self.find(self.parent[i]);
        }

        self.parent[i]
    }

    pub fn union(&mut self, x: usize, y: usize) {
        let rx = self.find(x);
        let ry = self.find(y);

        if rx != ry {
            if self.rank[rx] < self.rank[ry] {
                self.parent[rx] = ry;
            } else if self.rank[rx] > self.rank[ry] {
                self.parent[ry] = rx;
            } else {
                self.parent[ry] = rx;
                self.rank[rx] += 1;
            }
        }
    }
}

fn edit_distance(m1: &Motif, m2: &Motif) -> usize {
    // Line motifs according to modified base
    let mod_position_offset = m1.mod_position as i8 - m2.mod_position as i8;
    let mod_position_offset_abs = mod_position_offset.abs() as usize;
    let length_diff = m1.sequence.len() as i8 - m2.sequence.len() as i8;
    let length_diff_abs = length_diff.abs() as usize;

    let score = if mod_position_offset == 0 && length_diff == 0 {
        hamming_distance(&m1, &m2)
    } else if mod_position_offset != 0 && length_diff == 0 {
        // CCWG & CWGG
        let mut m1_extended = m1.clone();
        let mut m2_extended = m2.clone();

        if mod_position_offset > 0 {
            m1_extended.extend_motif_with_n(mod_position_offset_abs);
            m2_extended.prepend_n(mod_position_offset_abs);
        } else {
            m1_extended.prepend_n(mod_position_offset_abs);
            m2_extended.extend_motif_with_n(mod_position_offset_abs);
        }

        // N has a penalty of 0.5 but since an offset will always result in two Ns,
        // mod_position_offset_abs is just added.
        hamming_distance(&m1_extended, &m2_extended) + mod_position_offset_abs
    } else if mod_position_offset == 0 && length_diff != 0 {
        let mut m1_extended = m1.clone();
        let mut m2_extended = m2.clone();

        if length_diff > 0 {
            m2_extended.extend_motif_with_n(length_diff_abs);
        } else {
            m1_extended.extend_motif_with_n(length_diff_abs);
        }
        hamming_distance(&m1_extended, &m2_extended)
    } else {
        100
    };
    score
}

fn hamming_distance(s1: &Motif, s2: &Motif) -> usize {
    if s1.sequence.len() != s2.sequence.len() {
        panic!("Motif sequences should have the same length");
    }
    if s1.mod_position != s2.mod_position {
        panic!("Motifs should have the same mod_position");
    }
    if s1.mod_type != s2.mod_type {
        panic!("Motifs should have the same mod_type");
    }

    s1.sequence
        .iter()
        .zip(&s2.sequence)
        .fold(0, |score, (base1, base2)| {
            let possible_nucleotides_1 = base1.to_possible_nucleotides();
            let possible_nucleotides_2 = base2.to_possible_nucleotides();

            if possible_nucleotides_1
                .iter()
                .any(|x| possible_nucleotides_2.contains(x))
            {
                score
            } else {
                score + 1
            }
        })
}

fn cluster_motifs(motifs: &[Motif]) -> UnionFind {
    let n = motifs.len();
    let mut uf = UnionFind::new(n);

    for i in 0..n {
        for j in i + 1..n {
            if motifs[i].mod_type != motifs[j].mod_type {
                continue;
            } else if motifs[i].is_child_motif(&motifs[j])
                || motifs[j].is_child_motif(&motifs[i])
                || edit_distance(&motifs[i], &motifs[j]) <= 1
            {
                uf.union(i, j)
            }
        }
    }

    uf
}

fn group_motifs_by_set(uf: &mut UnionFind, motifs: &[Motif]) -> HashMap<usize, Vec<Motif>> {
    let n = motifs.len();
    let mut map: HashMap<usize, Vec<Motif>> = HashMap::new();

    for i in 0..n {
        let root = uf.find(i);
        map.entry(root).or_default().push(motifs[i].clone());
    }
    map
}

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

    let mut uf = cluster_motifs(&motifs);
    let motif_clusters = group_motifs_by_set(&mut uf, &motifs);

    // Within cluster find best candidate motif
    // Should be the smallest
    // In case several have the same lenght
    // they should be collapsed
    let mut motif_cluster_representatives = HashMap::new();

    for (_cluster, motifs_in_cluster) in motif_clusters {
        let min_motif = motifs_in_cluster
            .iter()
            .map(|m| m.sequence.len())
            .min()
            .unwrap();
        let smallest_motifs = motifs_in_cluster
            .iter()
            .cloned()
            .filter(|m| m.sequence.len() == min_motif)
            .collect::<Vec<_>>();

        let representative_motif = if smallest_motifs.len() > 1 {
            collapse_motifs(&smallest_motifs)?
        } else {
            smallest_motifs[0].clone()
        };

        motif_cluster_representatives.insert(representative_motif, motifs_in_cluster);
    }

    println!("{:#?}", motif_cluster_representatives);
    let outfile = std::fs::File::create(outpath)
        .with_context(|| format!("Failed to create file at: {:?}", outpath))?;
    let mut writer = BufWriter::new(outfile);

    writeln!(
        writer,
        "motif_representative\tmod_type_representative\tmod_position_representative\tmotif\tmod_type\tmod_position"
    )?;

    for (rep, motifs) in motif_cluster_representatives {
        let rep_motif_sequence = rep.sequence_to_string();
        let rep_mod_type_str = rep.mod_type.to_pileup_code();
        let rep_mod_position = rep.mod_position;

        for motif in motifs {
            let motif_sequence = motif.sequence_to_string();
            let mod_type_str = motif.mod_type.to_pileup_code();
            let mod_position = motif.mod_position;

            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                rep_motif_sequence,
                rep_mod_type_str,
                rep_mod_position,
                motif_sequence,
                mod_type_str,
                mod_position,
            )?;
            writer.flush()?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use methylome::Motif;

    use super::*;

    #[test]
    fn test_edit_distance_same_length_same_mod_pos() {
        let m1 = Motif::new("GATCC", "m", 3).unwrap();
        let m2 = Motif::new("GATCG", "m", 3).unwrap();
        let m3 = Motif::new("GTTCT", "m", 3).unwrap();

        let d1 = edit_distance(&m1, &m2);
        assert_eq!(d1, 1);
        let d2 = edit_distance(&m1, &m3);
        assert_eq!(d2, 2);
        let d3 = edit_distance(&m1, &m1);
        assert_eq!(d3, 0);
    }

    #[test]
    fn test_edit_distance_different_length_same_mod_pos() {
        let m1 = Motif::new("GATCC", "m", 3).unwrap();
        let m2 = Motif::new("GATC", "m", 3).unwrap();

        let d = edit_distance(&m1, &m2);
        assert_eq!(d, 0);
    }

    #[test]
    fn test_edit_distance_same_length_diff_mod_pos() {
        let m1 = Motif::new("CCWG", "m", 1).unwrap();
        let m2 = Motif::new("CWGG", "m", 0).unwrap();

        let d = edit_distance(&m1, &m2);
        assert_eq!(d, 1);
    }
    #[test]
    fn test_edit_distance_diff_length_diff_mod_pos() {
        let m1 = Motif::new("CCCCWG", "m", 1).unwrap();
        let m2 = Motif::new("CWGG", "m", 0).unwrap();

        let d = edit_distance(&m1, &m2);
        assert_eq!(d, 100);
    }

    // #[test]
    // fn test_levenshtein_distance_mod_position_unreachable_big_cost() {
    //     let m1 = Motif::new("CCWGG", "m", 1).unwrap();
    //     let m2 = Motif::new("GATCG", "m", 3).unwrap();

    //     let d = levenshtein_distance(&m1, &m2);
    //     assert_eq!(d, 3);
    // }
}
