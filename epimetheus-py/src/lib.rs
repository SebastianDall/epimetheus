use epimetheus_core::services::domain::parallel_processer::query_pileup;
use epimetheus_core::services::traits::PileupReader;
use methylome::{ModType, Strand};
use pyo3::prelude::*;
use std::path::Path;

use epimetheus_core::services::application::{
    methylation_pattern_service::extract_methylation_pattern,
    motif_clustering_service::motif_clustering,
};
use epimetheus_io::loaders::sequential_batch_loader::SequentialBatchLoader;
use epimetheus_io::readers::bedgz::Reader as GzPileupReader;
use epimetheus_io::readers::fasta::Reader as FastaReader;

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
            )?;

            meth_patthern.write_output(Path::new(output))
        })
    })
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
}

#[pyfunction]
fn remove_child_motifs(output: &str, motifs: Vec<String>) -> PyResult<()> {
    Python::with_gil(|py| py.allow_threads(|| motif_clustering(Path::new(output), &motifs)))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
}

// #[pyfunction]
// fn query_pileup_records(pileup_path: &str, contigs: Vec<String>) -> PyResult<Vec<String>> {
//     let mut reader = epimetheus_io::readers::bedgz::Reader::from_path(Path::new(pileup_path))
//         .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

//     let query_pileup_records = query_pileup(&mut reader, &contigs)
//         .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
// }

#[pymodule]
fn epymetheus(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(methylation_pattern, m)?)?;
    m.add_function(wrap_pyfunction!(remove_child_motifs, m)?)?;
    Ok(())
}
