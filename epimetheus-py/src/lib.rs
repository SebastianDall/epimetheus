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
use epimetheus_io::compression::bgzip::compressor::zip_pileup;
use pyo3::prelude::*;
use pyo3::types;
use std::path::Path;

use epimetheus_core::services::application::{
    methylation_pattern_service::extract_methylation_pattern,
    motif_clustering_service::motif_clustering,
};
use epimetheus_io::loaders::sequential_batch_loader::SequentialBatchLoader;
use epimetheus_io::readers::bgzf_bed::Reader as GzPileupReader;
use epimetheus_io::readers::fasta::Reader as FastaReader;

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

/// Query pileup records for specific contigs.
///
/// This function reads a pileup file and extracts all methylation records
/// for the specified contigs, returning them as a list of dictionaries.
///
/// Args:
///     pileup_path (str): Path to the pileup file (BED format, can be gzipped)
///     contigs (List[str]): List of contig names to query
///
/// Returns:
///     List[Dict]: List of dictionaries containing pileup record data with keys:
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
fn query_pileup_records(pileup_path: &str, contigs: Vec<String>) -> PyResult<PyObject> {
    let mut reader = epimetheus_io::readers::bgzf_bed::Reader::from_path(Path::new(pileup_path))
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let records = query_pileup(&mut reader, &contigs)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    Python::with_gil(|py| {
        let list = types::PyList::empty(py);

        for record in records {
            let dict = types::PyDict::new(py);
            dict.set_item("contig", record.contig)?;
            dict.set_item("start", record.start)?;
            dict.set_item("end", record.end)?;
            dict.set_item("mod_type", record.mod_type.to_pileup_code().to_string())?;
            dict.set_item("score", record.score)?;
            dict.set_item("strand", record.strand.to_string())?;
            dict.set_item("start_pos", record.start_pos)?;
            dict.set_item("end_pos", record.end_pos)?;
            dict.set_item("color", record.color)?;
            dict.set_item("n_valid_cov", record.n_valid_cov)?;
            dict.set_item("fraction_modified", record.fraction_modified)?;
            dict.set_item("n_modified", record.n_modified)?;
            dict.set_item("n_canonical", record.n_canonical)?;
            dict.set_item("n_other_mod", record.n_other_mod)?;
            dict.set_item("n_delete", record.n_delete)?;
            dict.set_item("n_fail", record.n_fail)?;
            dict.set_item("n_diff", record.n_diff)?;
            dict.set_item("n_no_call", record.n_no_call)?;

            list.append(dict)?;
        }
        Ok(list.into())
    })
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
    let output_path = output.map(Path::new);

    zip_pileup(Path::new(input), output_path, keep, force)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    Ok(())
}

#[pymodule]
fn epymetheus(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(methylation_pattern, m)?)?;
    m.add_function(wrap_pyfunction!(remove_child_motifs, m)?)?;
    m.add_function(wrap_pyfunction!(query_pileup_records, m)?)?;
    m.add_function(wrap_pyfunction!(bgzf_pileup, m)?)?;
    m.add_class::<MethylationOutput>()?;
    Ok(())
}
