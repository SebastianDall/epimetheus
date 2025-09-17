use criterion::{Criterion, black_box, criterion_group, criterion_main};
use std::{path::PathBuf, process::Command};

fn benchmark_methylation_pattern(c: &mut Criterion) {
    let mut group = c.benchmark_group("Methylation pattern");

    group.sample_size(10);

    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let data_dir = PathBuf::from(manifest_dir).join("tests/data");

    let pileup = data_dir.join("geobacillus.bed.gz");
    let assembly = data_dir.join("geobacillus-plasmids.assembly.fasta");

    // Test different batch sizes and thread counts
    let batch_sizes = [50, 100, 150];
    let thread_counts = [5, 10, 20];

    for &batch_size in &batch_sizes {
        for &threads in &thread_counts {
            let bench_name = format!(
                "methylation_pattern_batch_{}_threads_{}",
                batch_size, threads
            );

            group.bench_function(&bench_name, |b| {
                b.iter(|| {
                    let output_file = format!("target/bench_output_{}_{}.tsv", batch_size, threads);

                    let status = Command::new("cargo")
                        .args(&[
                            "run",
                            "--release",
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
                            &output_file,
                            "--batch-size",
                            &batch_size.to_string(),
                            "--threads",
                            &threads.to_string(),
                            "--min-valid-read-coverage",
                            "3",
                        ])
                        .status()
                        .expect("Failed to execute command");

                    black_box(status.success());

                    // Clean up output file
                    let _ = std::fs::remove_file(&output_file);
                });
            });
        }
    }
    group.finish();
}

criterion_group!(benches, benchmark_methylation_pattern);
criterion_main!(benches);
