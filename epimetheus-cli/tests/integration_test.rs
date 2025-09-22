use std::{fs, path::PathBuf, process::Command};
use tempfile::TempDir;

#[test]
fn test_methylation_pattern_median() {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let data_dir = PathBuf::from(manifest_dir).join("tests/data");

    let pileup = data_dir.join("geobacillus-plasmids.pileup.bed");
    let assembly = data_dir.join("geobacillus-plasmids.assembly.fasta");
    let expected_out = data_dir.join("expected_out_median.tsv");

    let out_file = PathBuf::from(manifest_dir)
        .join("target")
        .join("test_out_median.tsv");

    let status = Command::new("cargo")
        .args(&[
            "run",
            "--quiet",
            "--",
            "methylation-pattern",
            "-p",
            pileup.to_str().unwrap(),
            "-a",
            assembly.to_str().unwrap(),
            "-m",
            "GATC_a_1",
            "GATC_m_3",
            "RGATCY_a_2",
            "-o",
            out_file.to_str().unwrap(),
            "--batch-size",
            "2",
            "--min-valid-read-coverage",
            "3",
        ])
        .status()
        .expect("Failed to execute cargo run");

    assert!(
        status.success(),
        "Process ended with non-success status: {:?}",
        status
    );

    let actual = fs::read_to_string(&out_file).expect("Could not read output file");
    let expected = fs::read_to_string(&expected_out).expect("Could not read expected output file");

    let normalize = |s: &str| s.replace("\r\n", "\n");

    assert_eq!(
        normalize(actual.trim()),
        normalize(expected.trim()),
        "Output did not match expected"
    );
}

#[test]
fn test_methylation_pattern_weighted_mean() {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let data_dir = PathBuf::from(manifest_dir).join("tests/data");

    let pileup = data_dir.join("geobacillus-plasmids.pileup.bed");
    let assembly = data_dir.join("geobacillus-plasmids.assembly.fasta");
    let expected_out = data_dir.join("expected_out_weighted_mean.tsv");

    let out_file = PathBuf::from(manifest_dir)
        .join("target")
        .join("test_out_wm.tsv");

    let status = Command::new("cargo")
        .args(&[
            "run",
            "--quiet",
            "--",
            "methylation-pattern",
            "-p",
            pileup.to_str().unwrap(),
            "-a",
            assembly.to_str().unwrap(),
            "-m",
            "GATC_a_1",
            "GATC_m_3",
            "RGATCY_a_2",
            "-o",
            out_file.to_str().unwrap(),
            "--batch-size",
            "2",
            "--min-valid-read-coverage",
            "3",
            "--output-type",
            "weighted-mean",
        ])
        .status()
        .expect("Failed to execute cargo run");

    assert!(
        status.success(),
        "Process ended with non-success status: {:?}",
        status
    );

    let actual = fs::read_to_string(&out_file).expect("Could not read output file");
    let expected = fs::read_to_string(&expected_out).expect("Could not read expected output file");

    let normalize = |s: &str| s.replace("\r\n", "\n");

    assert_eq!(
        normalize(actual.trim()),
        normalize(expected.trim()),
        "Output did not match expected"
    );
}

#[test]
fn test_methylation_pattern_raw() {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let data_dir = PathBuf::from(manifest_dir).join("tests/data");

    let pileup = data_dir.join("geobacillus-plasmids.pileup.bed");
    let assembly = data_dir.join("geobacillus-plasmids.assembly.fasta");
    let expected_out = data_dir.join("expected_out_raw.tsv");

    let out_file = PathBuf::from(manifest_dir)
        .join("target")
        .join("test_out_raw.tsv");

    let status = Command::new("cargo")
        .args(&[
            "run",
            "--quiet",
            "--",
            "methylation-pattern",
            "-p",
            pileup.to_str().unwrap(),
            "-a",
            assembly.to_str().unwrap(),
            "-m",
            "GATC_a_1",
            "GATC_m_3",
            "RGATCY_a_2",
            "-o",
            out_file.to_str().unwrap(),
            "--batch-size",
            "2",
            "--min-valid-read-coverage",
            "3",
            "--output-type",
            "raw",
        ])
        .status()
        .expect("Failed to execute cargo run");

    assert!(
        status.success(),
        "Process ended with non-success status: {:?}",
        status
    );

    let actual = fs::read_to_string(&out_file).expect("Could not read output file");
    let expected = fs::read_to_string(&expected_out).expect("Could not read expected output file");

    let normalize = |s: &str| s.replace("\r\n", "\n");

    assert_eq!(
        normalize(actual.trim()),
        normalize(expected.trim()),
        "Output did not match expected"
    );
}

#[test]
fn test_methylation_pattern_median_gz() {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let data_dir = PathBuf::from(manifest_dir).join("tests/data");

    let pileup = data_dir.join("geobacillus.bed.gz");
    let assembly = data_dir.join("geobacillus-plasmids.assembly.fasta");
    let expected_out = data_dir.join("expected_out_median.tsv");

    let out_file = PathBuf::from(manifest_dir)
        .join("target")
        .join("test_out_median_gz.tsv");

    let status = Command::new("cargo")
        .args(&[
            "run",
            "--quiet",
            "--",
            "methylation-pattern",
            "-p",
            pileup.to_str().unwrap(),
            "-a",
            assembly.to_str().unwrap(),
            "-m",
            "GATC_a_1",
            "GATC_m_3",
            "RGATCY_a_2",
            "-o",
            out_file.to_str().unwrap(),
            "--batch-size",
            "2",
            "--min-valid-read-coverage",
            "3",
        ])
        .status()
        .expect("Failed to execute cargo run");

    assert!(
        status.success(),
        "Process ended with non-success status: {:?}",
        status
    );

    let actual = fs::read_to_string(&out_file).expect("Could not read output file");
    let expected = fs::read_to_string(&expected_out).expect("Could not read expected output file");

    let normalize = |s: &str| s.replace("\r\n", "\n");

    assert_eq!(
        normalize(actual.trim()),
        normalize(expected.trim()),
        "Output did not match expected"
    );
}

#[test]
fn test_methylation_pattern_weighted_mean_gz() {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let data_dir = PathBuf::from(manifest_dir).join("tests/data");

    let pileup = data_dir.join("geobacillus.bed.gz");
    let assembly = data_dir.join("geobacillus-plasmids.assembly.fasta");
    let expected_out = data_dir.join("expected_out_weighted_mean.tsv");

    let out_file = PathBuf::from(manifest_dir)
        .join("target")
        .join("test_out_wm_gz.tsv");

    let status = Command::new("cargo")
        .args(&[
            "run",
            "--quiet",
            "--",
            "methylation-pattern",
            "-p",
            pileup.to_str().unwrap(),
            "-a",
            assembly.to_str().unwrap(),
            "-m",
            "GATC_a_1",
            "GATC_m_3",
            "RGATCY_a_2",
            "-o",
            out_file.to_str().unwrap(),
            "--batch-size",
            "2",
            "--min-valid-read-coverage",
            "3",
            "--output-type",
            "weighted-mean",
        ])
        .status()
        .expect("Failed to execute cargo run");

    assert!(
        status.success(),
        "Process ended with non-success status: {:?}",
        status
    );

    let actual = fs::read_to_string(&out_file).expect("Could not read output file");
    let expected = fs::read_to_string(&expected_out).expect("Could not read expected output file");

    let normalize = |s: &str| s.replace("\r\n", "\n");

    assert_eq!(
        normalize(actual.trim()),
        normalize(expected.trim()),
        "Output did not match expected"
    );
}

#[test]
fn test_methylation_pattern_raw_gz() {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let data_dir = PathBuf::from(manifest_dir).join("tests/data");

    let pileup = data_dir.join("geobacillus.bed.gz");
    let assembly = data_dir.join("geobacillus-plasmids.assembly.fasta");
    let expected_out = data_dir.join("expected_out_raw.tsv");

    let out_file = PathBuf::from(manifest_dir)
        .join("target")
        .join("test_out_raw_gz.tsv");

    let status = Command::new("cargo")
        .args(&[
            "run",
            "--quiet",
            "--",
            "methylation-pattern",
            "-p",
            pileup.to_str().unwrap(),
            "-a",
            assembly.to_str().unwrap(),
            "-m",
            "GATC_a_1",
            "GATC_m_3",
            "RGATCY_a_2",
            "-o",
            out_file.to_str().unwrap(),
            "--batch-size",
            "2",
            "--min-valid-read-coverage",
            "3",
            "--output-type",
            "raw",
        ])
        .status()
        .expect("Failed to execute cargo run");

    assert!(
        status.success(),
        "Process ended with non-success status: {:?}",
        status
    );

    let actual = fs::read_to_string(&out_file).expect("Could not read output file");
    let expected = fs::read_to_string(&expected_out).expect("Could not read expected output file");

    let normalize = |s: &str| s.replace("\r\n", "\n");

    assert_eq!(
        normalize(actual.trim()),
        normalize(expected.trim()),
        "Output did not match expected"
    );
}

#[test]
fn test_compress_pileup() {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let data_dir = PathBuf::from(manifest_dir).join("tests/data");

    let pileup = data_dir.join("geobacillus-plasmids.pileup.bed");
    let temp_dir = TempDir::new().expect("Failed to create temp directory");
    let compressed_pileup = temp_dir.path().join("geobacillus-plasmids.pileup.bed.gz");

    let status = Command::new("cargo")
        .args(&[
            "run",
            "--quiet",
            "--",
            "bgzip",
            "compress",
            "-i",
            pileup.to_str().unwrap(),
            "-o",
            compressed_pileup.to_str().unwrap(),
            "--keep",
        ])
        .status()
        .expect("Failed to execute cargo run");

    assert!(
        status.success(),
        "Compression failed with status: {:?}",
        status
    );

    assert!(
        compressed_pileup.exists(),
        "Compressed file was not created: {:?}",
        compressed_pileup
    );

    assert!(
        temp_dir
            .path()
            .join("geobacillus-plasmids.pileup.bed.gz.tbi")
            .exists(),
        "Index file was not created"
    );
}

#[test]
fn test_verify_expected_outputs_from_raw() {
    use std::collections::HashMap;

    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let data_dir = PathBuf::from(manifest_dir).join("tests/data");

    let expected_raw = data_dir.join("expected_out_raw.tsv");
    let expected_weighted_mean = data_dir.join("expected_out_weighted_mean.tsv");
    let expected_median = data_dir.join("expected_out_median.tsv");

    // Read raw data
    let raw_data = fs::read_to_string(&expected_raw).expect("Could not read raw expected file");
    let raw_lines: Vec<&str> = raw_data.trim().lines().collect();

    // Parse header
    let header = raw_lines[0];
    assert_eq!(
        header,
        "contig\tstart\tstrand\tmotif\tmod_type\tmod_position\tn_modified\tn_valid_cov"
    );

    // Group data by (contig, motif, mod_type, mod_position)
    let mut grouped: HashMap<(String, String, String, u32), Vec<(u64, u64)>> = HashMap::new();

    for line in &raw_lines[1..] {
        let parts: Vec<&str> = line.split('\t').collect();
        let contig = parts[0].to_string();
        let motif = parts[3].to_string();
        let mod_type = parts[4].to_string();
        let mod_position: u32 = parts[5].parse().expect("Invalid mod_position");
        let n_modified: u64 = parts[6].parse().expect("Invalid n_modified");
        let n_valid_cov: u64 = parts[7].parse().expect("Invalid n_valid_cov");

        let key = (contig, motif, mod_type, mod_position);
        grouped
            .entry(key)
            .or_insert_with(Vec::new)
            .push((n_modified, n_valid_cov));
    }

    // Calculate weighted means and medians
    let mut calculated_weighted_means: HashMap<(String, String, String, u32), (f64, f64, u64)> =
        HashMap::new();
    let mut calculated_medians: HashMap<(String, String, String, u32), (f64, f64, u64)> =
        HashMap::new();

    for (key, values) in &grouped {
        let total_modified: u64 = values.iter().map(|(m, _)| m).sum();
        let total_coverage: u64 = values.iter().map(|(_, c)| c).sum();
        let motif_count = values.len() as u64;

        // Weighted mean calculation
        let weighted_mean_methylation = total_modified as f64 / total_coverage as f64;
        let mean_read_cov = total_coverage as f64 / motif_count as f64;

        calculated_weighted_means.insert(
            key.clone(),
            (weighted_mean_methylation, mean_read_cov, motif_count),
        );

        // Median calculation
        let mut methylation_ratios: Vec<f64> = values
            .iter()
            .map(|(modified, coverage)| *modified as f64 / *coverage as f64)
            .collect();
        methylation_ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let median_methylation = if methylation_ratios.len() % 2 == 0 {
            let mid = methylation_ratios.len() / 2;
            (methylation_ratios[mid - 1] + methylation_ratios[mid]) / 2.0
        } else {
            methylation_ratios[methylation_ratios.len() / 2]
        };

        calculated_medians.insert(
            key.clone(),
            (median_methylation, mean_read_cov, motif_count),
        );
    }

    // Read and verify weighted mean expected output
    let weighted_mean_data = fs::read_to_string(&expected_weighted_mean)
        .expect("Could not read weighted mean expected file");
    let weighted_mean_lines: Vec<&str> = weighted_mean_data.trim().lines().collect();

    for line in &weighted_mean_lines[1..] {
        let parts: Vec<&str> = line.split('\t').collect();
        let contig = parts[0].to_string();
        let motif = parts[1].to_string();
        let mod_type = parts[2].to_string();
        let mod_position: u32 = parts[3].parse().expect("Invalid mod_position");
        let expected_methylation: f64 = parts[4].parse().expect("Invalid methylation_value");
        let expected_mean_cov: f64 = parts[5].parse().expect("Invalid mean_read_cov");
        let expected_count: u64 = parts[6].parse().expect("Invalid n_motif_obs");

        let key = (contig, motif, mod_type, mod_position);
        let (calculated_methylation, calculated_mean_cov, calculated_count) =
            calculated_weighted_means.get(&key).expect(&format!(
                "Key {:?} not found in calculated weighted means",
                key
            ));

        assert!(
            (calculated_methylation - expected_methylation).abs() < 1e-10,
            "Weighted mean methylation mismatch for {:?}: calculated={}, expected={}",
            key,
            calculated_methylation,
            expected_methylation
        );

        assert!(
            (calculated_mean_cov - expected_mean_cov).abs() < 1e-10,
            "Mean coverage mismatch for {:?}: calculated={}, expected={}",
            key,
            calculated_mean_cov,
            expected_mean_cov
        );

        assert_eq!(
            *calculated_count, expected_count,
            "Motif count mismatch for {:?}: calculated={}, expected={}",
            key, calculated_count, expected_count
        );
    }

    // Read and verify median expected output
    let median_data =
        fs::read_to_string(&expected_median).expect("Could not read median expected file");
    let median_lines: Vec<&str> = median_data.trim().lines().collect();

    for line in &median_lines[1..] {
        let parts: Vec<&str> = line.split('\t').collect();
        let contig = parts[0].to_string();
        let motif = parts[1].to_string();
        let mod_type = parts[2].to_string();
        let mod_position: u32 = parts[3].parse().expect("Invalid mod_position");
        let expected_methylation: f64 = parts[4].parse().expect("Invalid methylation_value");
        let expected_mean_cov: f64 = parts[5].parse().expect("Invalid mean_read_cov");
        let expected_count: u64 = parts[6].parse().expect("Invalid n_motif_obs");

        let key = (contig, motif, mod_type, mod_position);
        let (calculated_methylation, calculated_mean_cov, calculated_count) = calculated_medians
            .get(&key)
            .expect(&format!("Key {:?} not found in calculated medians", key));

        assert!(
            (calculated_methylation - expected_methylation).abs() < 1e-10,
            "Median methylation mismatch for {:?}: calculated={}, expected={}",
            key,
            calculated_methylation,
            expected_methylation
        );

        assert!(
            (calculated_mean_cov - expected_mean_cov).abs() < 1e-10,
            "Mean coverage mismatch for {:?}: calculated={}, expected={}",
            key,
            calculated_mean_cov,
            expected_mean_cov
        );

        assert_eq!(
            *calculated_count, expected_count,
            "Motif count mismatch for {:?}: calculated={}, expected={}",
            key, calculated_count, expected_count
        );
    }
}
