//! Python bindings for Epimetheus methylation analysis toolkit.
//!
//! This module provides Python bindings for the core Epimetheus functionality,
//! allowing Python users to process methylation data from Nanopore sequencing.
//!
//! The main functions include:
//! - `methylation_pattern`: Extract methylation patterns for DNA motifs
//! - `remove_child_motifs`: Remove redundant child motifs through clustering
//! - `query_pileup_records`: Query specific contigs from pileup files
//! - `bgzf_pileup`: Compress pileup files using BGZF format

use epimetheus_core::models::methylation::MethylationOutput;
use epimetheus_core::services::domain::parallel_processer::query_pileup;
use epimetheus_core::services::traits::PileupReader;
use epimetheus_io::io::writers::bgzip::Writer;
use epimetheus_io::io::writers::bgzip::WriterType;
use epimetheus_io::services::compression_service::CompressorService;
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::path::PathBuf;

use epimetheus_core::services::application::{
    methylation_pattern_service::extract_methylation_pattern,
    motif_clustering_service::motif_clustering,
};
use epimetheus_io::io::readers::bed;
use epimetheus_io::io::readers::bgzf_bed::Reader as GzPileupReader;
use epimetheus_io::io::readers::fasta::Reader as FastaReader;
use epimetheus_io::loaders::sequential_batch_loader::SequentialBatchLoader;

/// Extract methylation patterns for specified DNA motifs from pileup data.
///
/// This function processes Nanopore methylation calls from a pileup file and extracts
/// methylation patterns for the specified DNA motifs, writing the results to an output file.
///
/// Args:
///     pileup (str): Path to the input pileup file (BED format, can be gzipped)
///     assembly (str): Path to the assembly FASTA file
///     output (str): Path for the output TSV file
///     threads (int): Number of threads to use for parallel processing
///     motifs (List[str]): List of DNA motifs to search for (e.g., ['GATC', 'CCWGG'])
///     min_valid_read_coverage (int): Minimum number of valid reads required for a position
///     batch_size (int): Number of records to process in each batch
///     min_valid_cov_to_diff_fraction (float): Minimum fraction of valid coverage to difference coverage
///     allow_assembly_pileup_mismatch (bool): Whether to allow mismatches between assembly and pileup
///     output_type (MethylationOutput): Output format type
///
/// Returns:
///     None
///
/// Raises:
///     PyRuntimeError: If processing fails due to IO errors or data format issues
#[pyfunction]
fn methylation_pattern(
    pileup: &str,
    assembly: &str,
    output: &str,
    threads: usize,
    motifs: Vec<String>,
    min_valid_read_coverage: usize,
    batch_size: usize,
    min_valid_cov_to_diff_fraction: f32,
    allow_assembly_pileup_mismatch: bool,
    output_type: MethylationOutput,
) -> PyResult<()> {
    Python::with_gil(|py| {
        py.allow_threads(|| {
            let meth_patthern = extract_methylation_pattern::<
                GzPileupReader,
                FastaReader,
                SequentialBatchLoader<std::io::BufReader<std::fs::File>>,
            >(
                Path::new(pileup),
                Path::new(assembly),
                // Path::new(output),
                threads,
                &motifs,
                min_valid_read_coverage as u32,
                batch_size,
                min_valid_cov_to_diff_fraction,
                allow_assembly_pileup_mismatch,
                &output_type,
            )?;

            meth_patthern.write_output(Path::new(output))
        })
    })
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
}

/// Remove child motifs from the output to avoid redundant patterns.
///
/// This function performs motif clustering to identify and remove child motifs
/// that are subsets of parent motifs, reducing redundancy in the results.
///
/// Args:
///     output (str): Path to the output file to process
///     motifs (List[str]): List of motifs to analyze for parent-child relationships
///
/// Returns:
///     None
///
/// Raises:
///     PyRuntimeError: If clustering fails due to IO errors or processing issues
#[pyfunction]
fn remove_child_motifs(output: &str, motifs: Vec<String>) -> PyResult<()> {
    Python::with_gil(|py| py.allow_threads(|| motif_clustering(Path::new(output), &motifs)))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
}

/// Query pileup records for specific contigs and return as Polars DataFrame.
///
/// This function reads a pileup file and extracts all methylation records
/// for the specified contigs, returning them as a Polars DataFrame for efficient processing.
///
/// Args:
///     pileup_path (str): Path to the pileup file (BED format, can be gzipped)
///     contigs (List[str]): List of contig names to query
///
/// Returns:
///     polars.DataFrame: DataFrame containing pileup record data with columns:
///         - contig: Contig/chromosome name
///         - start: Start position (0-based)
///         - end: End position
///         - mod_type: Modification type code
///         - score: Quality score
///         - strand: DNA strand (+ or -)
///         - start_pos: Start position in the feature
///         - end_pos: End position in the feature
///         - color: Color code
///         - n_valid_cov: Number of valid coverage reads
///         - fraction_modified: Fraction of reads showing modification
///         - n_modified: Number of modified reads
///         - n_canonical: Number of canonical reads
///         - n_other_mod: Number of other modification reads
///         - n_delete: Number of deletion reads
///         - n_fail: Number of failed reads
///         - n_diff: Number of different reads
///         - n_no_call: Number of no-call reads
///
/// Raises:
///     PyIOError: If the pileup file cannot be read
///     PyRuntimeError: If querying fails due to data processing issues
#[pyfunction]
fn query_pileup_records(pileup_path: &str, contigs: Vec<String>) -> PyResult<PyDataFrame> {
    let mut reader =
        epimetheus_io::io::readers::bgzf_bed::Reader::from_path(Path::new(pileup_path))
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    // Pre-allocate vectors for columns
    let mut contig_vec = Vec::new();
    let mut start_vec = Vec::new();
    let mut end_vec = Vec::new();
    let mut mod_type_vec = Vec::new();
    let mut score_vec = Vec::new();
    let mut strand_vec = Vec::new();
    let mut start_pos_vec = Vec::new();
    let mut end_pos_vec = Vec::new();
    let mut color_vec = Vec::new();
    let mut n_valid_cov_vec = Vec::new();
    let mut fraction_modified_vec = Vec::new();
    let mut n_modified_vec = Vec::new();
    let mut n_canonical_vec = Vec::new();
    let mut n_other_mod_vec = Vec::new();
    let mut n_delete_vec = Vec::new();
    let mut n_fail_vec = Vec::new();
    let mut n_diff_vec = Vec::new();
    let mut n_no_call_vec = Vec::new();

    for contig in contigs {
        let records = query_pileup(&mut reader, &[contig])
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        for record in records {
            contig_vec.push(record.contig);
            start_vec.push(record.start);
            end_vec.push(record.end);
            mod_type_vec.push(record.mod_type.to_pileup_code().to_string());
            score_vec.push(record.score);
            strand_vec.push(record.strand.to_string());
            start_pos_vec.push(record.start_pos);
            end_pos_vec.push(record.end_pos);
            color_vec.push(record.color);
            n_valid_cov_vec.push(record.n_valid_cov);
            fraction_modified_vec.push(record.fraction_modified);
            n_modified_vec.push(record.n_modified);
            n_canonical_vec.push(record.n_canonical);
            n_other_mod_vec.push(record.n_other_mod);
            n_delete_vec.push(record.n_delete);
            n_fail_vec.push(record.n_fail);
            n_diff_vec.push(record.n_diff);
            n_no_call_vec.push(record.n_no_call);
        }
    }

    // Create DataFrame from columns
    let df = df! [
        "contig" => contig_vec,
        "start" => start_vec,
        "end" => end_vec,
        "mod_type" => mod_type_vec,
        "score" => score_vec,
        "strand" => strand_vec,
        "start_pos" => start_pos_vec,
        "end_pos" => end_pos_vec,
        "color" => color_vec,
        "n_valid_cov" => n_valid_cov_vec,
        "fraction_modified" => fraction_modified_vec,
        "n_modified" => n_modified_vec,
        "n_canonical" => n_canonical_vec,
        "n_other_mod" => n_other_mod_vec,
        "n_delete" => n_delete_vec,
        "n_fail" => n_fail_vec,
        "n_diff" => n_diff_vec,
        "n_no_call" => n_no_call_vec,
    ]
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    Ok(PyDataFrame(df))
}

/// Compress a pileup file using BGZF compression.
///
/// This function compresses a pileup file using the BGZF (Blocked GZip Format)
/// compression algorithm, which allows for efficient random access.
///
/// Args:
///     input (str): Path to the input pileup file
///     output (Optional[str]): Path for the compressed output file.
///                           If None, adds .gz extension to input filename
///     keep (bool): Whether to keep the original uncompressed file
///     force (bool): Whether to overwrite existing output file
///
/// Returns:
///     None
///
/// Raises:
///     PyRuntimeError: If compression fails due to IO errors or file access issues
#[pyfunction]
fn bgzf_pileup(input: &str, output: Option<&str>, keep: bool, force: bool) -> PyResult<()> {
    let input_path = Path::new(input);
    let input_file = File::open(input_path)
        .map_err(|e| pyo3::exceptions::PyFileNotFoundError::new_err(e.to_string()))?;

    let output_str = match output {
        Some(out) => out,
        None => &format!("{}.gz", input.to_string()),
    };
    let output_path = Path::new(output_str);

    if !keep {
        Python::with_gil(|py| {
            let warnings = py.import("warnings")?;
            let message = format!(
                "File removal is enabled, will remove '{}' after compression. Set 'keep' to remove this behavior.",
                &input_path.display()
            );

            warnings.call_method1("warn", (message,))?;
            Ok::<(), PyErr>(())
        })?;
    }
    if !force & output_path.exists() {
        let message = format!(
            "File '{}' already exists. Set '--force' to override.",
            &output_path.display()
        );

        // warnings.call_method1("warn", (message,))?;
        return Err(pyo3::exceptions::PyFileExistsError::new_err(message));
    }

    let reader = bed::InputReader::File(bed::LineReader::new(BufReader::new(input_file)));

    CompressorService::compress_pileup(reader, Some(output_path))
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    if !keep {
        Python::with_gil(|py| {
            let warnings = py.import("warnings")?;
            let message = format!("Removing original file: {}.", &input_path.display());

            warnings.call_method1("warn", (message,))?;

            std::fs::remove_file(&input_path)?;
            Ok::<(), PyErr>(())
        })?;
    }
    Ok(())
}

#[pyclass]
pub struct BgzfWriter {
    writer: Option<WriterType>,
    output_path: PathBuf,
}

#[pymethods]
impl BgzfWriter {
    #[new]
    fn new(output_path: &str, force: bool) -> PyResult<Self> {
        let path = PathBuf::from(output_path);

        if !force && path.exists() {
            return Err(pyo3::exceptions::PyFileExistsError::new_err(format!(
                "File: '{}' already exists",
                output_path
            )));
        }

        let writer = WriterType::File(
            Writer::from_path(&path)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?,
        );

        Ok(Self {
            writer: Some(writer),
            output_path: path,
        })
    }

    fn write_lines(&mut self, lines: Vec<String>) -> PyResult<()> {
        if let Some(ref mut writer) = self.writer {
            writer
                .compress_from_lines(lines.into_iter())
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        }

        Ok(())
    }

    fn finish(&mut self) -> PyResult<()> {
        if let Some(mut writer) = self.writer.take() {
            let tbx_path = format!("{}.tbi", self.output_path.display());
            writer
                .write_tabix(Path::new(&tbx_path))
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

            writer
                .finish()
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        }
        Ok(())
    }
}

#[pymodule]
fn epymetheus(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(methylation_pattern, m)?)?;
    m.add_function(wrap_pyfunction!(remove_child_motifs, m)?)?;
    m.add_function(wrap_pyfunction!(query_pileup_records, m)?)?;
    m.add_function(wrap_pyfunction!(bgzf_pileup, m)?)?;
    // m.add_function(wrap_pyfunction!(bgzf_pileup_from_lines, m)?)?;
    m.add_class::<MethylationOutput>()?;
    m.add_class::<BgzfWriter>()?;
    Ok(())
}
