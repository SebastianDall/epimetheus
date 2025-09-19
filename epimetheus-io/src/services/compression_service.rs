use anyhow::Result;
use std::path::Path;

use crate::io::{
    readers::bed::InputReader,
    writers::bgzip::{Writer, WriterType},
};

pub struct CompressorService;

impl CompressorService {
    pub fn compress_pileup(input_reader: InputReader, output: Option<&Path>) -> Result<()> {
        let mut writer = match output {
            Some(path) => WriterType::File(Writer::from_path(path)?),
            None => WriterType::StdOut(Writer::to_stdout()?),
        };

        match input_reader {
            InputReader::File(reader) => writer.compress_from_reader(reader)?,
            InputReader::StdIn(reader) => writer.compress_from_reader(reader)?,
        }

        if let Some(path) = output {
            let tbx_path = format!("{}.tbi", path.display());
            writer.write_tabix(Path::new(&tbx_path))?;
        }

        writer.finish()?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_zip_pileup_creates_gz_and_tbi_files() {
        // Create a temporary input file with BED data
        let mut input_file = NamedTempFile::new().unwrap();
        writeln!(input_file, "contig_1\t0\t100\tA\t50").unwrap();
        writeln!(input_file, "contig_1\t150\t250\tC\t75").unwrap();
        writeln!(input_file, "contig_2\t300\t400\tG\t25").unwrap();
        input_file.flush().unwrap();

        // Create temporary output directory
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("test_output.bed.gz");

        // Test the compression function
        let result = zip_pileup(input_file.path(), Some(output_path.as_path()), true, false);
        assert!(result.is_ok(), "zip_pileup failed: {:?}", result.err());

        // Check that the compressed file was created
        assert!(output_path.exists(), "Output .gz file was not created");

        // Check that the tabix index was created
        let tbi_path = format!("{}.tbi", output_path.display());
        assert!(
            Path::new(&tbi_path).exists(),
            "Tabix .tbi file was not created"
        );

        // Verify the compressed file can be read
        let compressed_reader = File::open(&output_path).map(bgzf::io::Reader::new).unwrap();
        let mut buf_reader = BufReader::new(compressed_reader);
        let mut line = String::new();
        let mut line_count = 0;

        while buf_reader.read_line(&mut line).unwrap() > 0 {
            line_count += 1;
            line.clear();
        }

        assert_eq!(line_count, 3, "Compressed file should contain 3 lines");

        // Verify the original file still exists (because keep=true)
        assert!(input_file.path().exists(), "Original file should be kept");
    }

    #[test]
    fn test_zip_pileup_removes_original_when_keep_false() {
        // Create a real temporary file (not NamedTempFile which auto-deletes)
        let temp_dir = tempfile::tempdir().unwrap();
        let input_path = temp_dir.path().join("input.bed");
        let mut input_file = File::create(&input_path).unwrap();
        writeln!(input_file, "chr1\t0\t100\tA\t50").unwrap();
        input_file.sync_all().unwrap();
        drop(input_file);

        let output_path = temp_dir.path().join("output.bed.gz");

        let result = zip_pileup(&input_path, Some(output_path.as_path()), false, false);
        assert!(result.is_ok(), "zip_pileup failed: {:?}", result.err());

        // Check that the original file was removed
        assert!(
            !input_path.exists(),
            "Original file should have been removed"
        );

        // Check that output files were created
        assert!(output_path.exists(), "Output .gz file was not created");
        let tbi_path = format!("{}.tbi", output_path.display());
        assert!(
            Path::new(&tbi_path).exists(),
            "Tabix .tbi file was not created"
        );
    }

    #[test]
    fn test_zip_pileup_handles_zero_coordinates() {
        let mut input_file = NamedTempFile::new().unwrap();
        writeln!(input_file, "chr1\t0\t1\tA\t50").unwrap(); // Zero start coordinate
        input_file.flush().unwrap();

        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("zero_coord.bed.gz");

        let result = zip_pileup(input_file.path(), Some(output_path.as_path()), true, false);
        assert!(
            result.is_ok(),
            "zip_pileup should handle zero coordinates: {:?}",
            result.err()
        );

        assert!(output_path.exists(), "Output file should be created");
        let tbi_path = format!("{}.tbi", output_path.display());
        assert!(
            Path::new(&tbi_path).exists(),
            "Tabix index should be created"
        );
    }

    #[test]
    fn test_zip_pileup_validates_output_extension() {
        let mut input_file = NamedTempFile::new().unwrap();
        writeln!(input_file, "chr1\t0\t100\tA\t50").unwrap();
        input_file.flush().unwrap();

        let result = zip_pileup(
            input_file.path(),
            Some(Path::new("invalid_extension.txt")),
            true,
            false,
        );
        assert!(result.is_err(), "Should fail with invalid extension");

        let error_msg = result.unwrap_err().to_string();
        assert!(
            error_msg.contains("must have .gz extension"),
            "Error should mention .gz extension requirement"
        );
    }
}
