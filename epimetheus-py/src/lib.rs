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

#[pyfunction]
fn remove_child_motifs(output: &str, motifs: Vec<String>) -> PyResult<()> {
    Python::with_gil(|py| py.allow_threads(|| motif_clustering(Path::new(output), &motifs)))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
}

#[pyfunction]
fn query_pileup_records(pileup_path: &str, contigs: Vec<String>) -> PyResult<PyObject> {
    let mut reader = epimetheus_io::readers::bedgz::Reader::from_path(Path::new(pileup_path))
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
    Ok(())
}
