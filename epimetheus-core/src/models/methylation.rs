use anyhow::{Result, bail};
use methylome::{ModType, Motif, Strand};

#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct MethylationCoverage {
    n_modified: u32,
    n_valid_cov: u32,
}

impl MethylationCoverage {
    pub fn new(n_modified: u32, n_valid_cov: u32) -> Result<Self> {
        if n_modified > n_valid_cov {
            bail!(
                "Invalid coverage: n_valid_cov ({}) cannot be less than n_modified ({})",
                n_valid_cov,
                n_modified
            )
        }

        Ok(Self {
            n_modified,
            n_valid_cov,
        })
    }

    // pub fn get_n_modified(&self) -> u32 {
    //     self.n_modified
    // }

    pub fn get_n_valid_cov(&self) -> u32 {
        self.n_valid_cov
    }

    pub fn fraction_modified(&self) -> f64 {
        self.n_modified as f64 / self.n_valid_cov as f64
    }
}

pub struct MethylationRecord {
    pub contig: String,
    pub position: usize,
    pub strand: Strand,
    pub mod_type: ModType,
    pub methylation: MethylationCoverage,
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

    #[allow(dead_code)]
    pub fn get_contig_id(&self) -> String {
        self.contig.to_string()
    }
}

pub struct MotifMethylationDegree {
    pub contig: String,
    pub motif: Motif,
    pub median: f64,
    pub mean_read_cov: f64,
    pub n_motif_obs: u32,
    pub motif_occurences_total: u32,
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_methylation_coverage_valid() -> Result<()> {
        // Test valid inputs
        let coverage = MethylationCoverage::new(5, 10)?;
        assert_eq!(coverage.n_modified, 5);
        assert_eq!(coverage.n_valid_cov, 10);

        let coverage = MethylationCoverage::new(0, 0)?;
        assert_eq!(coverage.n_modified, 0);
        assert_eq!(coverage.n_valid_cov, 0);

        Ok(())
    }

    #[test]
    fn test_methylation_coverage_invalid() {
        // Test invalid input: n_valid_cov < n_modified
        let result = MethylationCoverage::new(10, 5);

        assert!(result.is_err());
        if let Err(e) = result {
            assert_eq!(
                e.to_string(),
                "Invalid coverage: n_valid_cov (5) cannot be less than n_modified (10)"
            );
        }
    }
}
