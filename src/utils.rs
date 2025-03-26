use anyhow::{anyhow, Context, Result};
use std::{fs, path::Path};

pub fn create_output_file(outpath: &Path) -> Result<()> {
    if let Some(ext) = outpath.extension() {
        if ext != "tsv" {
            anyhow::bail!("Incorrect file extension {:?}. Should be tsv", ext);
        }
        Ok(if let Some(parent) = outpath.parent() {
            fs::create_dir_all(parent)
                .with_context(|| format!("Could not create parent directory: {:?}", parent))?;
        })
    } else {
        return Err(anyhow!(
            "No filename provided for output. Should be a .tsv file."
        ));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_create_output_file() {
        let dir = tempdir().unwrap();
        let file = dir.path().join("output.tsv");
        let file_path = Path::new(file.to_str().unwrap());

        let result = create_output_file(file_path);

        assert!(result.is_ok(), "File creation should succeed");

        assert!(
            file.parent().unwrap().exists(),
            "Parent directory should exist"
        );

        assert!(dir.path().exists(), "Temporary directory should exist");
    }

    #[test]
    fn test_create_output_file_incorrect_extension() {
        let dir = tempdir().unwrap();
        let file = dir.path().join("output.csv");
        let file_path = Path::new(file.to_str().unwrap());

        let result = create_output_file(file_path);

        assert!(result.is_err(), "File creation should fail");
        assert_eq!(
            result.unwrap_err().to_string(),
            "Incorrect file extension \"csv\". Should be tsv"
        );
    }
}
