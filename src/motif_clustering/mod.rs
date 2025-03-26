use anyhow::{Context, Result};
use log::info;
use methylome::Motif;
use std::{collections::HashMap, path::Path};

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

fn levenshtein_distance(s1: &Motif, s2: &Motif) -> usize {
    let m = s1.sequence.len();
    let n = s2.sequence.len();

    // let mod1 = s1.mod_position as usize;
    // let mod2 = s2.mod_position as usize;

    let mut dp: Vec<Vec<usize>> = vec![vec![0; n + 1]; m + 1];
    // let big_cost = 9999 as usize;

    // Populate the matrix sides
    for i in 0..(m + 1) {
        dp[i][0] = i
    }
    for j in 0..(n + 1) {
        dp[0][j] = j
    }

    for i in 1..(m + 1) {
        for j in 1..(n + 1) {
            let p1 = s1.sequence[i - 1].to_possible_nucleotides();
            let p2 = s2.sequence[j - 1].to_possible_nucleotides();
            let cost = {
                if p1.iter().any(|x| p2.contains(x)) {
                    0
                } else {
                    1
                }
            };

            let del = dp[i - 1][j] + 1;
            let ins = dp[i][j - 1] + 1;
            let sub = dp[i - 1][j - 1] + cost;

            // if (i - 1 == mod1) && (j - 1 != mod2) {
            //     del = big_cost;
            //     ins = big_cost;
            //     sub = big_cost;
            // } else if (j - 1 == mod2) && (i - 1 != mod1) {
            //     del = big_cost;
            //     ins = big_cost;
            //     sub = big_cost;
            // }

            dp[i][j] = *[del, ins, sub].iter().min().unwrap();
        }
    }
    dp[m][n]
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
                || levenshtein_distance(&motifs[i], &motifs[j]) <= 1
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

pub fn motif_clustering(args: MotifClusteringArgs) -> Result<()> {
    let output_path = Path::new(&args.output);

    create_output_file(output_path)?;

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
    let motif_cluster_representatives = HashMap::new();

    for (cluster, motifs_in_cluster) in motif_clusters {
        let min_motif = motifs_in_cluster
            .iter()
            .map(|m| m.sequence.len())
            .min()
            .unwrap();
        let smallest_motifs = motifs_in_cluster
            .iter()
            .filter(|m| m.sequence.len() == min_motif)
            .collect::<Vec<_>>();

        println!("{}: {:#?}", cluster, smallest_motifs);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use methylome::Motif;

    use super::*;

    #[test]
    fn test_levenshtein_distance_same_length() {
        let m1 = Motif::new("GATCC", "m", 3).unwrap();
        let m2 = Motif::new("GATCG", "m", 3).unwrap();
        let m3 = Motif::new("GTTCT", "m", 3).unwrap();

        let d1 = levenshtein_distance(&m1, &m2);
        assert_eq!(d1, 1);
        let d2 = levenshtein_distance(&m1, &m3);
        assert_eq!(d2, 2);
        let d3 = levenshtein_distance(&m1, &m1);
        assert_eq!(d3, 0);
    }

    #[test]
    fn test_levenshtein_distance_different_length() {
        let m1 = Motif::new("GATCC", "m", 3).unwrap();
        let m2 = Motif::new("GATC", "m", 3).unwrap();

        let d = levenshtein_distance(&m1, &m2);
        assert_eq!(d, 0);
    }

    #[test]
    fn test_levenshtein_distance_ambiguous_bases() {
        let m1 = Motif::new("GATCR", "m", 3).unwrap();
        let m2 = Motif::new("GATCG", "m", 3).unwrap();

        let d = levenshtein_distance(&m1, &m2);
        assert_eq!(d, 0);
    }

    // #[test]
    // fn test_levenshtein_distance_mod_position_unreachable_big_cost() {
    //     let m1 = Motif::new("CCWGG", "m", 1).unwrap();
    //     let m2 = Motif::new("GATCG", "m", 3).unwrap();

    //     let d = levenshtein_distance(&m1, &m2);
    //     assert_eq!(d, 3);
    // }
}
