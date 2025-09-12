use anyhow::{Result, anyhow};
use csv::StringRecord;
use methylome::{ModType, Strand};

use crate::models::methylation::{MethylationCoverage, MethylationRecord};

pub fn parse_to_methylation_record(
    contig: String,
    record: &StringRecord,
    min_valid_read_coverage: u32,
    min_valid_cov_to_diff_fraction: f32,
) -> Result<Option<MethylationRecord>> {
    let position: usize = record
        .get(1)
        .ok_or_else(|| anyhow!("Missing position field."))?
        .parse()
        .map_err(|_| anyhow!("Invalid position field"))?;

    let mod_type: ModType = record
        .get(3)
        .ok_or_else(|| anyhow!("Missing modification type field."))?
        .parse()?;

    let strand: Strand = record
        .get(5)
        .ok_or_else(|| anyhow!("Missing strand field"))?
        .parse()?;

    let n_valid_cov: u32 = record
        .get(9)
        .ok_or_else(|| anyhow!("Missing n_valid_cov field."))?
        .parse()
        .map_err(|_| anyhow!("Invalid n_valid_cov field."))?;

    if n_valid_cov < min_valid_read_coverage {
        return Ok(None);
    }
    let n_modified: u32 = record
        .get(11)
        .ok_or_else(|| anyhow!("Missing n_modified field."))?
        .parse()
        .map_err(|_| anyhow!("Invalid n_modified field"))?;

    let n_other_mod: u32 = record
        .get(13)
        .ok_or_else(|| anyhow!("Missing n_other_mod field."))?
        .parse()
        .map_err(|_| anyhow!("Invalid n_other_mod field."))?;

    if n_other_mod > n_modified {
        return Ok(None);
    }

    let n_diff: u32 = record
        .get(16)
        .ok_or_else(|| anyhow!("Missing n_diff field."))?
        .parse()
        .map_err(|_| anyhow!("Invalid n_diff field."))?;

    // n_diff is the number of reads with another base than the canonical. The fraction of valid cov
    // compared to total cov should be higher thatn the threshold.
    if (n_valid_cov as f32 / (n_diff as f32 + n_valid_cov as f32)) < min_valid_cov_to_diff_fraction
    {
        return Ok(None);
    }

    let methylation = MethylationCoverage::new(n_modified, n_valid_cov - n_other_mod)?;

    let methylation_record =
        MethylationRecord::new(contig, position, strand, mod_type, methylation);

    Ok(Some(methylation_record))
}
