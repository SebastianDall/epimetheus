use pyo3::prelude::*;



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
    allow_assembly_pilup_mismatch: bool,
) -> PyResult<()> {
    let args = epimetheus_core::extract_methylation_pattern::MethylationPatternArgs {
        pileup: pileup.to_string(),
        assembly: assembly.to_string(),
        output: output.to_string(),
        threads,
        motifs: Some(motifs),
        min_valid_read_coverage: min_valid_read_coverage as u32,
        batch_size,
        min_valid_cov_to_diff_fraction,
        allow_assembly_pilup_mismatch,
    };


    
    Python::with_gil(|py| {
        py.allow_threads(|| {
            epimetheus_core::extract_methylation_pattern(args)
        })
    })
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
}

#[pyfunction]
fn remove_child_motifs(
        output: &str,
        motifs: Vec<String>,
) -> PyResult<()> {
    let args = epimetheus_core::motif_clustering::MotifClusteringArgs {
        output: output.to_string(),
        motifs: Some(motifs),
    };


    
    Python::with_gil(|py| {
        py.allow_threads(|| {
            epimetheus_core::motif_clustering(args)
        })
    })
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
}

#[pymodule]
fn epymetheus(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(methylation_pattern, m)?)?;
    m.add_function(wrap_pyfunction!(remove_child_motifs, m)?)?;
    Ok(())
}
