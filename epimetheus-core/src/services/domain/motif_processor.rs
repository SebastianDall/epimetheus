use methylome::Motif;
use anyhow::Context;
use std::str::FromStr;


pub fn create_motifs(motifs_str: &Vec<String>) -> anyhow::Result<Vec<Motif>> {
    motifs_str.into_iter().map(|motif| {
        let parts: Vec<&str> = motif.split("_").collect();

        if parts.len() != 3 {
            anyhow::bail!(
                "Invalid motif format '{}' encountered. Expected format: '<sequence>_<mod_type>_<mod_position>'",
                motif
            );
        }

            let sequence = parts[0];
            let mod_type = parts[1];
            let mod_position = u8::from_str(parts[2]).with_context(|| {
                format!("Failed to parse mod_position '{}' in motif '{}'.", parts[2], motif)
            })?;

            Motif::new(sequence, mod_type, mod_position).with_context(|| {
                format!("Failed to create motif from '{}'", motif)
            })
        
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_motifs_success() {
        let motifs_args = vec!["GATC_a_1".to_string()];
        let result = create_motifs(&motifs_args);
        assert!(
            result.is_ok(),
            "Expected Ok, but got err: {:?}",
            result.err()
        );
    }
    #[test]
    fn test_create_motifs_failure() {
        let motifs_args = vec!["GATC_a_3".to_string()];
        let result = create_motifs(&motifs_args);
        assert!(
            result.is_err(),
            "Expected Err, but got Ok: {:?}",
            result.ok()
        );
    }
}
