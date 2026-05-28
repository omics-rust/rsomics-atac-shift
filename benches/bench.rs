use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_atac_shift(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-atac-shift");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bam = manifest.join("tests/golden/bench_small.bam");
    let out = tempfile::NamedTempFile::new().unwrap();

    c.bench_function("rsomics-atac-shift golden", |b| {
        b.iter(|| {
            let status = Command::new(black_box(bin))
                .args([bam.to_str().unwrap(), "-o", out.path().to_str().unwrap()])
                .status()
                .unwrap();
            assert!(status.success());
        });
    });
}

criterion_group!(benches, bench_atac_shift);
criterion_main!(benches);
