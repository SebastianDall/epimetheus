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
            InputReader::Lines(lines) => writer.compress_from_lines(lines)?,
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
    use crate::io::readers::bed::{InputReader, LineReader};
    use noodles_bgzf as bgzf;
    use std::{
        fs::File,
        io::{BufRead, BufReader, Write},
    };
    use tempfile::NamedTempFile;

    fn create_test_bed_data() -> NamedTempFile {
        let mut input_file = NamedTempFile::new().unwrap();
        writeln!(
            input_file,
            "contig_3\t0\t1\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )
        .unwrap();
        writeln!(
            input_file,
            "contig_3\t6\t7\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )
        .unwrap();
        writeln!(
            input_file,
            "contig_3\t10\t11\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )
        .unwrap();
        input_file.flush().unwrap();
        input_file
    }

    #[test]
    fn test_compress_pileup_creates_gz_and_tbi_files() {
        let input_file = create_test_bed_data();
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("test_output.bed.gz");

        // Create InputReader from file
        let file = File::open(input_file.path()).unwrap();
        let line_reader = LineReader::new(BufReader::new(file));
        let input_reader = InputReader::File(line_reader);

        // Test compression
        let result = CompressorService::compress_pileup(input_reader, Some(&output_path));
        assert!(result.is_ok(), "compress_pileup failed: {:?}", result.err());

        // Verify outputs
        assert!(output_path.exists(), "Output .gz file was not created");
        let tbi_path = format!("{}.tbi", output_path.display());
        assert!(
            Path::new(&tbi_path).exists(),
            "Tabix .tbi file was not created"
        );

        // Verify compressed content
        let compressed_reader = File::open(&output_path).map(bgzf::io::Reader::new).unwrap();
        let mut buf_reader = BufReader::new(compressed_reader);
        let mut line = String::new();
        let mut line_count = 0;

        while buf_reader.read_line(&mut line).unwrap() > 0 {
            line_count += 1;
            line.clear();
        }
        assert_eq!(line_count, 3, "Compressed file should contain 3 lines");
    }

    #[test]
    fn test_compress_pileup_to_stdout() {
        let input_file = create_test_bed_data();

        // Create InputReader from file
        let file = File::open(input_file.path()).unwrap();
        let line_reader = LineReader::new(BufReader::new(file));
        let input_reader = InputReader::File(line_reader);

        // Test compression to stdout (no output path)
        let result = CompressorService::compress_pileup(input_reader, None);
        assert!(
            result.is_ok(),
            "compress_pileup to stdout failed: {:?}",
            result.err()
        );
        // Note: Can't easily test stdout output, but we verify it doesn't crash
    }

    #[test]
    fn test_compress_handles_zero_coordinates() {
        let mut input_file = NamedTempFile::new().unwrap();
        writeln!(
            input_file,
            "contig_3\t0\t1\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )
        .unwrap(); // Zero start coordinate
        input_file.flush().unwrap();

        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("zero_coord.bed.gz");

        let file = File::open(input_file.path()).unwrap();
        let line_reader = LineReader::new(BufReader::new(file));
        let input_reader = InputReader::File(line_reader);

        let result = CompressorService::compress_pileup(input_reader, Some(&output_path));
        assert!(
            result.is_ok(),
            "Should handle zero coordinates: {:?}",
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
    fn test_compress_from_memory_data() {
        // Create a temporary file with test data instead of using Cursor
        let mut input_file = NamedTempFile::new().unwrap();
        writeln!(
            input_file,
            "contig_3\t0\t1\ta\t133\t+\t0\t1\t255,0,0\t15\t0.00\t15\t123\t0\t0\t6\t0\t0"
        )
        .unwrap();
        input_file.flush().unwrap();

        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("memory_test.bed.gz");

        let file = File::open(input_file.path()).unwrap();
        let line_reader = LineReader::new(BufReader::new(file));
        let input_reader = InputReader::File(line_reader);

        let result = CompressorService::compress_pileup(input_reader, Some(&output_path));
        assert!(
            result.is_ok(),
            "Should handle file input: {:?}",
            result.err()
        );

        assert!(output_path.exists(), "Output file should be created");
    }
}
