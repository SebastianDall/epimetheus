use std::{path::Path, str::FromStr};

use ahash::AHashMap;
use anyhow::{Result, bail};
use clap::ValueEnum;
use methylome::{ModType, Motif, Strand};

#[cfg(feature = "python")]
use pyo3::{IntoPyObject, types::PyAnyMethods};

use crate::models::{
    contig::{ContigId, Position as ContigPosition},
    pileup::PileupRecord,
};

#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct MethylationCoverage {
    n_modified: u32,
    n_valid_cov: u32,
}

impl MethylationCoverage {
    pub fn new(n_modified: u32, n_valid_cov: u32, n_other_mod: u32) -> Result<Self> {
        if n_modified > n_valid_cov {
            bail!(
                "Invalid coverage: n_valid_cov ({}) cannot be less than n_modified ({})",
                n_valid_cov,
                n_modified
            )
        }

        let n_valid_cov = n_valid_cov - n_other_mod;

        Ok(Self {
            n_modified,
            n_valid_cov,
        })
    }

    pub fn get_n_modified(&self) -> u32 {
        self.n_modified
    }

    pub fn get_n_valid_cov(&self) -> u32 {
        self.n_valid_cov
    }

    pub fn fraction_modified(&self) -> f64 {
        self.n_modified as f64 / self.n_valid_cov as f64
    }
}

#[derive(Clone)]
pub struct MethylationRecord {
    pub contig: String,
    pub position: usize,
    pub strand: Strand,
    pub mod_type: ModType,
    pub methylation: MethylationCoverage,
}

impl TryFrom<PileupRecord> for MethylationRecord {
    type Error = anyhow::Error;

    fn try_from(value: PileupRecord) -> std::result::Result<Self, Self::Error> {
        let meth =
            MethylationCoverage::new(value.n_modified, value.n_valid_cov, value.n_other_mod)?;
        Ok(Self {
            contig: value.contig.clone(),
            position: value.start as usize,
            strand: value.strand.clone(),
            mod_type: value.mod_type.clone(),
            methylation: meth,
        })
    }
}
impl MethylationRecord {
    pub fn new(
        contig: String,
        position: usize,
        strand: Strand,
        mod_type: ModType,
        methylation: MethylationCoverage,
    ) -> Self {
        Self {
            contig,
            position,
            strand,
            mod_type,
            methylation,
        }
    }

    pub fn try_from_with_filters(
        value: PileupRecord,
        min_valid_read_coverage: u32,
        min_valid_cov_to_diff_fraction: f32,
    ) -> Result<Option<Self>> {
        if value.n_valid_cov < min_valid_read_coverage {
            return Ok(None);
        }

        if value.n_other_mod > value.n_modified {
            return Ok(None);
        }

        if (value.n_valid_cov as f32 / (value.n_diff as f32 + value.n_valid_cov as f32))
            < min_valid_cov_to_diff_fraction
        {
            return Ok(None);
        }

        let meth =
            MethylationCoverage::new(value.n_modified, value.n_valid_cov, value.n_other_mod)?;

        Ok(Some(Self {
            contig: value.contig.clone(),
            position: value.start as usize,
            strand: value.strand.clone(),
            mod_type: value.mod_type.clone(),
            methylation: meth,
        }))
    }

    #[allow(dead_code)]
    pub fn get_contig_id(&self) -> String {
        self.contig.to_string()
    }
}

pub trait MotifMethylationDegree {
    fn get_contig(&self) -> &str;
    fn get_motif(&self) -> &Motif;
    fn get_methylation_value(&self) -> f64;
    fn get_mean_read_cov(&self) -> f64;
    fn get_n_motif_obs(&self) -> u32;
    // fn get_motif_occurences_total(&self) -> u32;

    fn to_csv_line(&self, delim: char) -> String {
        let motif_seq = self.get_motif().sequence_to_string();
        let mod_type = self.get_motif().mod_type.to_pileup_code();
        let mod_position = self.get_motif().mod_position;

        format!(
            "{}{delim}{}{delim}{}{delim}{}{delim}{}{delim}{}{delim}{}",
            self.get_contig(),
            motif_seq,
            mod_type,
            mod_position,
            self.get_methylation_value(),
            self.get_mean_read_cov(),
            self.get_n_motif_obs(),
            // self.get_motif_occurences_total(),
        )
    }
}

#[derive(PartialEq, Clone, PartialOrd)]
pub struct MedianMotifMethylationDegree {
    pub contig: String,
    pub motif: Motif,
    pub median: f64,
    pub mean_read_cov: f64,
    pub n_motif_obs: u32,
    // pub motif_occurences_total: u32,
}

impl MotifMethylationDegree for MedianMotifMethylationDegree {
    fn get_contig(&self) -> &str {
        self.contig.as_str()
    }

    fn get_motif(&self) -> &Motif {
        &self.motif
    }

    fn get_methylation_value(&self) -> f64 {
        self.median
    }

    fn get_mean_read_cov(&self) -> f64 {
        self.mean_read_cov
    }

    fn get_n_motif_obs(&self) -> u32 {
        self.n_motif_obs
    }

    // fn get_motif_occurences_total(&self) -> u32 {
    //     self.motif_occurences_total
    // }
}

#[derive(PartialEq, Clone, PartialOrd)]
pub struct WeightedMeanMotifMethylationDegree {
    pub contig: String,
    pub motif: Motif,
    pub w_mean: f64,
    pub mean_read_cov: f64,
    pub n_motif_obs: u32,
    // pub motif_occurences_total: u32,
}

impl MotifMethylationDegree for WeightedMeanMotifMethylationDegree {
    fn get_contig(&self) -> &str {
        self.contig.as_str()
    }

    fn get_motif(&self) -> &Motif {
        &self.motif
    }

    fn get_methylation_value(&self) -> f64 {
        self.w_mean
    }

    fn get_mean_read_cov(&self) -> f64 {
        self.mean_read_cov
    }

    fn get_n_motif_obs(&self) -> u32 {
        self.n_motif_obs
    }

    // fn get_motif_occurences_total(&self) -> u32 {
    //     self.motif_occurences_total
    // }
}

pub struct MotifMethylationPositions {
    pub methylation: AHashMap<(ContigId, Motif, ContigPosition, Strand), MethylationCoverage>,
}

impl MotifMethylationPositions {
    pub fn new(
        methylation: AHashMap<(ContigId, Motif, ContigPosition, Strand), MethylationCoverage>,
    ) -> Self {
        Self { methylation }
    }

    fn group_by_motif(&self) -> AHashMap<(ContigId, Motif), Vec<&MethylationCoverage>> {
        let mut grouped: AHashMap<(ContigId, Motif), Vec<&MethylationCoverage>> = AHashMap::new();

        for ((contig_id, motif, _position, _strand), coverage) in &self.methylation {
            grouped
                .entry((contig_id.clone(), motif.clone()))
                .or_insert_with(Vec::new)
                .push(coverage);
        }

        grouped
    }

    pub fn to_median_degrees(&self) -> Vec<MedianMotifMethylationDegree> {
        self.group_by_motif()
            .into_iter()
            .map(|((contig_id, motif), coverages)| {
                let mut fractions: Vec<f64> = coverages
                    .iter()
                    .map(|cov| cov.fraction_modified())
                    .collect();

                fractions.sort_by(|a, b| a.partial_cmp(b).unwrap());

                let median = if fractions.len() % 2 == 0 {
                    let mid = fractions.len() / 2;
                    (fractions[mid - 1] + fractions[mid]) / 2.0
                } else {
                    fractions[fractions.len() / 2]
                };

                let mean_read_cov = {
                    let total_cov: u64 = coverages
                        .iter()
                        .map(|cov| cov.get_n_valid_cov() as u64)
                        .sum();

                    total_cov as f64 / coverages.len() as f64
                };

                MedianMotifMethylationDegree {
                    contig: contig_id,
                    motif,
                    median,
                    mean_read_cov,
                    n_motif_obs: coverages.len() as u32,
                }
            })
            .collect()
    }

    pub fn to_weighted_mean_degress(&self) -> Vec<WeightedMeanMotifMethylationDegree> {
        self.group_by_motif()
            .into_iter()
            .map(|((contig_id, motif), coverages)| {
                let fraction_weight = coverages
                    .iter()
                    .map(|cov| cov.fraction_modified() * cov.get_n_valid_cov() as f64)
                    .sum::<f64>();

                let total_weights = coverages
                    .iter()
                    .map(|cov| cov.get_n_valid_cov())
                    .sum::<u32>();

                let weighted_mean = fraction_weight / total_weights as f64;

                let mean_read_cov = {
                    let total_cov: u64 = coverages
                        .iter()
                        .map(|cov| cov.get_n_valid_cov() as u64)
                        .sum();

                    total_cov as f64 / coverages.len() as f64
                };

                WeightedMeanMotifMethylationDegree {
                    contig: contig_id,
                    motif,
                    w_mean: weighted_mean,
                    mean_read_cov,
                    n_motif_obs: coverages.len() as u32,
                }
            })
            .collect()
    }
}

#[derive(Debug, Clone, ValueEnum)]
#[cfg_attr(feature = "python", pyo3::pyclass(module = "epymetheus"))]
pub enum MethylationOutput {
    Raw,
    Median,
    WeightedMean,
}

impl ToString for MethylationOutput {
    fn to_string(&self) -> String {
        match self {
            Self::Raw => "raw".to_string(),
            Self::Median => "median".to_string(),
            Self::WeightedMean => "weighted_mean".to_string(),
        }
    }
}

impl FromStr for MethylationOutput {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "raw" => Ok(Self::Raw),
            "median" => Ok(Self::Median),
            "weighted_mean" => Ok(Self::WeightedMean),
            _ => Err(format!("Invalid output type: {}", s)),
        }
    }
}

#[cfg(feature = "python")]
#[pyo3::pymethods]
impl MethylationOutput {
    fn __reduce__(&self) -> pyo3::PyResult<(pyo3::PyObject, (String,))> {
        pyo3::Python::with_gil(|py| {
            let state = match self {
                MethylationOutput::Median => "Median".to_string(),
                MethylationOutput::WeightedMean => "WeightedMean".to_string(),
                MethylationOutput::Raw => "Raw".to_string(),
            };
            let constructor = py.get_type::<MethylationOutput>().getattr("_from_state")?;
            Ok((constructor.into_pyobject(py)?.unbind(), (state,)))
        })
    }

    #[classmethod]
    fn _from_state(
        _cls: &pyo3::Bound<'_, pyo3::types::PyType>,
        state: String,
    ) -> pyo3::PyResult<MethylationOutput> {
        match state.as_str() {
            "Median" => Ok(MethylationOutput::Median),
            "WeightedMean" => Ok(MethylationOutput::WeightedMean),
            "Raw" => Ok(MethylationOutput::Raw),
            _ => Err(pyo3::PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Invalid state",
            )),
        }
    }

    fn __eq__(&self, other: &Self) -> bool {
        std::mem::discriminant(self) == std::mem::discriminant(other)
    }
}

pub enum MethylationPatternVariant {
    Raw(MotifMethylationPositions),
    Median(Vec<MedianMotifMethylationDegree>),
    WeightedMean(Vec<WeightedMeanMotifMethylationDegree>),
}

impl MethylationPatternVariant {
    pub fn write_output<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        use std::fs::File;
        use std::io::{BufWriter, Write};

        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        match self {
            MethylationPatternVariant::Raw(meth_pos) => {
                writeln!(
                    writer,
                    "contig\tstart\tstrand\tmotif\tmod_type\tmod_position\tn_modified\tn_valid_cov"
                )?;

                let mut sorted_entries: Vec<_> = meth_pos.methylation.iter().collect();
                sorted_entries.sort_by_key(|((contig_id, motif, pos, strand), _)| {
                    (contig_id.clone(), motif.clone(), *pos, strand)
                });

                for ((contig_id, motif, pos, strand), meth) in sorted_entries {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        contig_id,
                        pos,
                        strand.to_string(),
                        motif.sequence_to_string(),
                        motif.mod_type.to_pileup_code(),
                        motif.mod_position,
                        meth.get_n_modified(),
                        meth.get_n_valid_cov(),
                    )?;
                }
            }
            MethylationPatternVariant::Median(degrees) => {
                writeln!(
                    writer,
                    "contig\tmotif\tmod_type\tmod_position\tmethylation_value\tmean_read_cov\tn_motif_obs"
                )?;
                let mut sorted_degrees = degrees.clone();
                sorted_degrees.sort_by(|a, b| a.partial_cmp(b).expect("Ordering failed"));

                for deg in sorted_degrees {
                    writeln!(writer, "{}", deg.to_csv_line('\t'))?;
                }
            }
            MethylationPatternVariant::WeightedMean(degrees) => {
                writeln!(
                    writer,
                    "contig\tmotif\tmod_type\tmod_position\tmethylation_value\tmean_read_cov\tn_motif_obs"
                )?;
                let mut sorted_degrees = degrees.clone();
                sorted_degrees.sort_by(|a, b| a.partial_cmp(b).expect("Ordering failed"));

                for deg in sorted_degrees {
                    writeln!(writer, "{}", deg.to_csv_line('\t'))?;
                }
            }
        }

        writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_methylation_coverage_valid() -> Result<()> {
        // Test valid inputs
        let coverage = MethylationCoverage::new(5, 10, 0)?;
        assert_eq!(coverage.n_modified, 5);
        assert_eq!(coverage.n_valid_cov, 10);

        let coverage = MethylationCoverage::new(0, 0, 0)?;
        assert_eq!(coverage.n_modified, 0);
        assert_eq!(coverage.n_valid_cov, 0);

        Ok(())
    }

    #[test]
    fn test_methylation_coverage_invalid() {
        // Test invalid input: n_valid_cov < n_modified
        let result = MethylationCoverage::new(10, 5, 0);

        assert!(result.is_err());
        if let Err(e) = result {
            assert_eq!(
                e.to_string(),
                "Invalid coverage: n_valid_cov (5) cannot be less than n_modified (10)"
            );
        }
    }
}
