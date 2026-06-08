//! Compatibility against `deeptools alignmentSieve --ATACshift` v3.5.6.
//!
//! Two flavours. The `*_matches_golden` tests diff ours against a committed
//! golden captured once from alignmentSieve; they decode both BAMs in-process
//! with noodles and always run in CI (no samtools, no live oracle). The
//! live-oracle tests re-shift through alignmentSieve itself and self-skip when
//! deeptools isn't installed.
//!
//! Comparison is field-level (name, flag, pos, cigar, tlen, pnext), sorted by
//! (name, flag): byte-exact BAM diff is impossible since deeptools stamps a
//! different @PG header and may emit records in a different order.
use std::io::Write as _;
use std::path::PathBuf;
use std::process::Command;

use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::io::Write as _;

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

fn alignmentsieve_available() -> bool {
    Command::new("alignmentSieve")
        .arg("--version")
        .output()
        .is_ok()
}

/// Parse headerless SAM lines into (name, flag, pos_1based, cigar, tlen, next_pos_1based) tuples,
/// sorted by (name, flag) so output-order differences between tools don't produce false mismatches.
fn parse_sam_fields_sorted(output: &str) -> Vec<(String, u16, i64, String, i32, i64)> {
    let mut records: Vec<(String, u16, i64, String, i32, i64)> = output
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| {
            let fields: Vec<&str> = l.split('\t').collect();
            let name = fields[0].to_string();
            let flag: u16 = fields[1].parse().unwrap();
            let pos: i64 = fields[3].parse().unwrap(); // 1-based SAM
            let cigar = fields[5].to_string();
            let tlen: i32 = fields[8].parse().unwrap();
            let next_pos: i64 = fields[7].parse().unwrap(); // 1-based
            (name, flag, pos, cigar, tlen, next_pos)
        })
        .collect();
    records.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    records
}

/// Decode a BAM into headerless SAM record lines, the same shape `samtools view`
/// emits, so the records can be field-diffed without spawning samtools.
fn bam_to_sam_text(bam: &std::path::Path) -> String {
    let mut reader = bam::io::reader::Builder.build_from_path(bam).unwrap();
    let header = reader.read_header().unwrap();

    let mut out = Vec::new();
    let mut writer = sam::io::Writer::new(&mut out);
    for result in reader.records() {
        writer
            .write_alignment_record(&header, &result.unwrap())
            .unwrap();
    }
    writer.get_mut().flush().unwrap();
    String::from_utf8(out).unwrap()
}

type SamFields = (String, u16, i64, String, i32, i64);

fn run_ours(input_bam: &std::path::Path, out: &std::path::Path) {
    let binary = env!("CARGO_BIN_EXE_rsomics-atac-shift");
    let status = Command::new(binary)
        .args([
            input_bam.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "-q",
        ])
        .status()
        .expect("rsomics-atac-shift failed to launch");
    assert!(status.success(), "rsomics-atac-shift exited non-zero");
}

fn assert_fields_eq(ours: &[SamFields], expect: &[SamFields], oracle: &str) {
    assert_eq!(
        ours.len(),
        expect.len(),
        "record count mismatch: ours={} {oracle}={}",
        ours.len(),
        expect.len()
    );
    for (i, (a, b)) in ours.iter().zip(expect.iter()).enumerate() {
        assert_eq!(a.0, b.0, "record {i}: name {:?} vs {:?}", a.0, b.0);
        assert_eq!(a.1, b.1, "record {i} ({}): flag {} vs {}", a.0, a.1, b.1);
        assert_eq!(a.2, b.2, "record {i} ({}): POS {} vs {}", a.0, a.2, b.2);
        assert_eq!(
            a.3, b.3,
            "record {i} ({}): CIGAR {:?} vs {:?}",
            a.0, a.3, b.3
        );
        assert_eq!(a.4, b.4, "record {i} ({}): TLEN {} vs {}", a.0, a.4, b.4);
        assert_eq!(a.5, b.5, "record {i} ({}): PNEXT {} vs {}", a.0, a.5, b.5);
    }
}

fn run_compat_check(input_bam: &std::path::Path) {
    let tmp = tempfile::tempdir().expect("tempdir");
    let our_out = tmp.path().join("our_shifted.bam");
    let dt_out = tmp.path().join("dt_shifted.bam");

    run_ours(input_bam, &our_out);

    let dt_status = Command::new("alignmentSieve")
        .args([
            "-b",
            input_bam.to_str().unwrap(),
            "-o",
            dt_out.to_str().unwrap(),
            "--ATACshift",
            "--numberOfProcessors",
            "1",
        ])
        .status()
        .expect("alignmentSieve failed to launch");
    assert!(dt_status.success(), "alignmentSieve exited non-zero");

    let our_records = parse_sam_fields_sorted(&bam_to_sam_text(&our_out));
    let dt_records = parse_sam_fields_sorted(&bam_to_sam_text(&dt_out));
    assert_fields_eq(&our_records, &dt_records, "deeptools");
}

/// Diff ours against a committed golden produced once by `alignmentSieve
/// --ATACshift` v3.5.6. Always runs in CI: both BAMs are decoded in-process via
/// noodles, so neither deeptools nor samtools need be installed.
fn run_golden_check(input_bam: &std::path::Path, golden_bam: &std::path::Path) {
    let tmp = tempfile::tempdir().expect("tempdir");
    let our_out = tmp.path().join("our_shifted.bam");
    run_ours(input_bam, &our_out);

    let our_records = parse_sam_fields_sorted(&bam_to_sam_text(&our_out));
    let golden_records = parse_sam_fields_sorted(&bam_to_sam_text(golden_bam));
    assert_fields_eq(&our_records, &golden_records, "golden");
}

#[test]
fn atac_shift_matches_deeptools() {
    if !alignmentsieve_available() {
        eprintln!("SKIP: alignmentSieve not found (deeptools not installed)");
        return;
    }

    let fixtures = fixture_dir();
    let input_bam = fixtures.join("atac_small.bam");
    if !input_bam.exists() {
        eprintln!("SKIP: fixture BAM not found at {:?}", input_bam);
        return;
    }

    run_compat_check(&input_bam);
}

/// Representative compat test using a ~3700-record fixture derived from bench.bam
/// (chr1:1-30000, ~65% soft-clipped reads). This exercises the CIGAR class that
/// the 50M-only atac_small.bam fixture cannot catch: reads with leading or trailing
/// soft clips (e.g. `3S97M`) where deeptools uses `query_alignment_end` (query
/// coordinate) rather than the reference-consuming CIGAR span to compute `end`,
/// producing a different new-CIGAR length (96M vs 93M for `3S97M` with +4 shift).
#[test]
fn atac_shift_softclip_compat() {
    if !alignmentsieve_available() {
        eprintln!("SKIP: alignmentSieve not found (deeptools not installed)");
        return;
    }

    let fixtures = fixture_dir();
    let input_bam = fixtures.join("bench_small.bam");
    if !input_bam.exists() {
        eprintln!("SKIP: bench_small.bam not found");
        return;
    }

    run_compat_check(&input_bam);
}

#[test]
fn atac_shift_matches_golden() {
    let fixtures = fixture_dir();
    run_golden_check(
        &fixtures.join("atac_small.bam"),
        &fixtures.join("golden_shifted.bam"),
    );
}

#[test]
fn atac_shift_softclip_matches_golden() {
    let fixtures = fixture_dir();
    run_golden_check(
        &fixtures.join("bench_small.bam"),
        &fixtures.join("golden_bench_shifted.bam"),
    );
}
