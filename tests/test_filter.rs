use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use tempfile::TempDir;

/// A simple FASTQ record for building test data.
struct FqRecord {
    name: &'static str,
    seq: &'static str,
    qual: &'static str,
}

/// Writes a kraken2 report file and returns its path.
fn write_report(dir: &Path, lines: &[&str]) -> PathBuf {
    let path = dir.join("report.txt");
    let content = lines.join("\n") + "\n";
    std::fs::write(&path, content).unwrap();
    path
}

/// Writes an empty (0-byte) kraken2 report file and returns its path.
fn write_empty_report(dir: &Path) -> PathBuf {
    let path = dir.join("report.txt");
    std::fs::write(&path, "").unwrap();
    path
}

/// Writes a kraken2 per-read output file and returns its path.
fn write_kraken_output(dir: &Path, records: &[&str]) -> PathBuf {
    let path = dir.join("kraken_output.txt");
    let content = records.join("\n") + "\n";
    std::fs::write(&path, content).unwrap();
    path
}

/// Writes a FASTQ file (plain text) and returns its path.
fn write_fastq(dir: &Path, name: &str, records: &[FqRecord]) -> PathBuf {
    let path = dir.join(name);
    let mut f = File::create(&path).unwrap();
    for rec in records {
        writeln!(f, "@{}\n{}\n+\n{}", rec.name, rec.seq, rec.qual).unwrap();
    }
    path
}

/// Writes a gzip-compressed FASTQ file and returns its path.
fn write_fastq_gz(dir: &Path, name: &str, records: &[FqRecord]) -> PathBuf {
    let path = dir.join(name);
    let file = File::create(&path).unwrap();
    let mut encoder = GzEncoder::new(file, Compression::fast());
    for rec in records {
        writeln!(encoder, "@{}\n{}\n+\n{}", rec.name, rec.seq, rec.qual).unwrap();
    }
    encoder.finish().unwrap();
    path
}

/// Parses FASTQ lines and returns (name, seq) pairs.
fn parse_fastq_lines(lines: &[String]) -> Vec<(String, String)> {
    let mut result = Vec::new();
    for chunk in lines.chunks(4) {
        assert!(chunk[0].starts_with('@'), "expected FASTQ header, got: {}", chunk[0]);
        let name = chunk[0][1..].to_string();
        let seq = chunk[1].clone();
        result.push((name, seq));
    }
    result
}

/// Reads a bgzf/gzip FASTQ file and returns (name, seq) pairs.
fn read_output_fastq(path: &Path) -> Vec<(String, String)> {
    let file = File::open(path).unwrap();
    let reader = BufReader::new(MultiGzDecoder::new(file));
    let lines: Vec<String> = reader.lines().map(|l| l.unwrap()).collect();
    parse_fastq_lines(&lines)
}

/// Parses uncompressed FASTQ bytes (e.g. captured stdout) into (name, seq) pairs.
fn parse_plain_fastq(bytes: &[u8]) -> Vec<(String, String)> {
    let text = std::str::from_utf8(bytes).unwrap();
    let lines: Vec<String> = text.lines().map(str::to_string).collect();
    parse_fastq_lines(&lines)
}

/// Returns true if the file starts with the gzip magic bytes.
fn is_gzip_file(path: &Path) -> bool {
    std::fs::read(path).unwrap().starts_with(&[0x1f, 0x8b])
}

/// Returns the path to the k2tools binary.
fn k2tools_bin() -> PathBuf {
    let mut path = std::env::current_exe().unwrap();
    path.pop(); // remove test binary name
    path.pop(); // remove 'deps'
    path.push("k2tools");
    path
}

/// Runs k2tools filter with the given arguments and returns the output.
fn run_filter(args: &[&str]) -> std::process::Output {
    Command::new(k2tools_bin()).arg("filter").args(args).output().unwrap()
}

// A standard report with:
// unclassified(0): 100 reads
// root(1): 900 reads, 5 direct
//   Bacteria(2): 600 reads, 100 direct
//     E.coli(3): 500 reads, 500 direct
//   Eukaryota(4): 300 reads, 100 direct
//     Human(5): 200 reads, 200 direct
fn standard_report_lines() -> Vec<&'static str> {
    vec![
        " 10.00\t100\t100\tU\t0\tunclassified",
        " 90.00\t900\t5\tR\t1\troot",
        " 60.00\t600\t100\tD\t2\t  Bacteria",
        " 50.00\t500\t500\tS\t3\t    Escherichia coli",
        " 30.00\t300\t100\tD\t4\t  Eukaryota",
        " 20.00\t200\t200\tS\t5\t    Homo sapiens",
    ]
}

fn make_reads() -> Vec<FqRecord> {
    vec![
        FqRecord { name: "r1", seq: "ACGT", qual: "IIII" }, // taxon 3 (E.coli)
        FqRecord { name: "r2", seq: "TGCA", qual: "IIII" }, // unclassified
        FqRecord { name: "r3", seq: "AAAA", qual: "IIII" }, // taxon 5 (Human)
        FqRecord { name: "r4", seq: "CCCC", qual: "IIII" }, // taxon 2 (Bacteria direct)
        FqRecord { name: "r5", seq: "GGGG", qual: "IIII" }, // taxon 3 (E.coli)
    ]
}

/// R2 mates for `make_reads()`, with sequences distinct from their R1 mates.
fn make_reads_r2() -> Vec<FqRecord> {
    vec![
        FqRecord { name: "r1", seq: "TTTT", qual: "IIII" },
        FqRecord { name: "r2", seq: "AAAA", qual: "IIII" },
        FqRecord { name: "r3", seq: "CCCC", qual: "IIII" },
        FqRecord { name: "r4", seq: "GGGG", qual: "IIII" },
        FqRecord { name: "r5", seq: "TTTT", qual: "IIII" },
    ]
}

fn make_kraken_lines() -> Vec<&'static str> {
    vec![
        "C\tr1\t3\t4\t3:1",
        "U\tr2\t0\t4\t0:1",
        "C\tr3\t5\t4\t5:1",
        "C\tr4\t2\t4\t2:1",
        "C\tr5\t3\t4\t3:1",
    ]
}

#[test]
fn test_single_end_exact_taxon() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let input = write_fastq(dir.path(), "input.fq", &make_reads());
    let output = dir.path().join("output.fq.gz");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "3",
        "--threads",
        "2",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    let records = read_output_fastq(&output);
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].0, "r1");
    assert_eq!(records[0].1, "ACGT");
    assert_eq!(records[1].0, "r5");
    assert_eq!(records[1].1, "GGGG");
}

#[test]
fn test_paired_end() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());

    let in1 = write_fastq(dir.path(), "r1.fq", &make_reads());
    let in2 = write_fastq(dir.path(), "r2.fq", &make_reads_r2());
    let out1 = dir.path().join("out1.fq.gz");
    let out2 = dir.path().join("out2.fq.gz");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        in1.to_str().unwrap(),
        in2.to_str().unwrap(),
        "-o",
        out1.to_str().unwrap(),
        out2.to_str().unwrap(),
        "-t",
        "5",
        "--threads",
        "2",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));

    let recs1 = read_output_fastq(&out1);
    let recs2 = read_output_fastq(&out2);
    assert_eq!(recs1.len(), 1);
    assert_eq!(recs2.len(), 1);
    assert_eq!(recs1[0].0, "r3");
    assert_eq!(recs1[0].1, "AAAA");
    assert_eq!(recs2[0].0, "r3");
    assert_eq!(recs2[0].1, "CCCC");
}

#[test]
fn test_with_descendants() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let input = write_fastq(dir.path(), "input.fq", &make_reads());
    let output = dir.path().join("output.fq.gz");

    // Bacteria (2) + descendants (E.coli 3)
    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "2",
        "-d",
        "--threads",
        "2",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    let records = read_output_fastq(&output);
    // r1 (taxon 3), r4 (taxon 2), r5 (taxon 3)
    assert_eq!(records.len(), 3);
    assert_eq!(records[0].0, "r1");
    assert_eq!(records[1].0, "r4");
    assert_eq!(records[2].0, "r5");
}

#[test]
fn test_include_unclassified_only() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let input = write_fastq(dir.path(), "input.fq", &make_reads());
    let output = dir.path().join("output.fq.gz");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-u",
        "--threads",
        "2",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    let records = read_output_fastq(&output);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].0, "r2");
}

#[test]
fn test_taxon_ids_plus_unclassified() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let input = write_fastq(dir.path(), "input.fq", &make_reads());
    let output = dir.path().join("output.fq.gz");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "5",
        "-u",
        "--threads",
        "2",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    let records = read_output_fastq(&output);
    // r2 (unclassified) + r3 (Human taxon 5)
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].0, "r2");
    assert_eq!(records[1].0, "r3");
}

#[test]
fn test_no_matches() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let input = write_fastq(dir.path(), "input.fq", &make_reads());
    let output = dir.path().join("output.fq.gz");

    // taxon 4 (Eukaryota) has no direct assignments in the kraken output
    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "4",
        "--threads",
        "2",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    let records = read_output_fastq(&output);
    assert_eq!(records.len(), 0);
}

#[test]
fn test_all_match() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let input = write_fastq(dir.path(), "input.fq", &make_reads());
    let output = dir.path().join("output.fq.gz");

    // Root with descendants covers everything except unclassified; add -u too
    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "1",
        "-d",
        "-u",
        "--threads",
        "2",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    let records = read_output_fastq(&output);
    assert_eq!(records.len(), 5);
}

#[test]
fn test_multiple_taxon_ids() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let input = write_fastq(dir.path(), "input.fq", &make_reads());
    let output = dir.path().join("output.fq.gz");

    // taxon 3 (E.coli) + taxon 5 (Human)
    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "3",
        "5",
        "--threads",
        "2",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    let records = read_output_fastq(&output);
    // r1 (3), r3 (5), r5 (3)
    assert_eq!(records.len(), 3);
    assert_eq!(records[0].0, "r1");
    assert_eq!(records[1].0, "r3");
    assert_eq!(records[2].0, "r5");
}

#[test]
fn test_gzip_compressed_input() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let input = write_fastq_gz(dir.path(), "input.fq.gz", &make_reads());
    let output = dir.path().join("output.fq.gz");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "3",
        "--threads",
        "2",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    let records = read_output_fastq(&output);
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].0, "r1");
    assert_eq!(records[1].0, "r5");
}

#[test]
fn test_fastq_longer_than_kraken() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());

    // Only 3 kraken lines but 5 FASTQ records
    let kraken = write_kraken_output(
        dir.path(),
        &["C\tr1\t3\t4\t3:1", "U\tr2\t0\t4\t0:1", "C\tr3\t5\t4\t5:1"],
    );
    let input = write_fastq(dir.path(), "input.fq", &make_reads());
    let output = dir.path().join("output.fq.gz");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "3",
        "--threads",
        "2",
    ]);

    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("more records than kraken"), "unexpected error: {stderr}");
}

#[test]
fn test_kraken_longer_than_fastq() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());

    // Only 3 FASTQ records but 5 kraken lines
    let reads = vec![
        FqRecord { name: "r1", seq: "ACGT", qual: "IIII" },
        FqRecord { name: "r2", seq: "TGCA", qual: "IIII" },
        FqRecord { name: "r3", seq: "AAAA", qual: "IIII" },
    ];
    let input = write_fastq(dir.path(), "input.fq", &reads);
    let output = dir.path().join("output.fq.gz");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "3",
        "--threads",
        "2",
    ]);

    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("FASTQ input ended before kraken"), "unexpected error: {stderr}");
}

#[test]
fn test_mismatched_input_output_count() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let input = write_fastq(dir.path(), "input.fq", &make_reads());
    let out1 = dir.path().join("out1.fq.gz");
    let out2 = dir.path().join("out2.fq.gz");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        out1.to_str().unwrap(),
        out2.to_str().unwrap(),
        "-t",
        "3",
    ]);

    assert!(!result.status.success());
}

#[test]
fn test_mismatched_read_names() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());

    // Kraken says read names are r1, r2 but FASTQ has x1, x2
    let kraken = write_kraken_output(dir.path(), &["C\tr1\t3\t4\t3:1", "C\tr2\t3\t4\t3:1"]);
    let reads = vec![
        FqRecord { name: "x1", seq: "ACGT", qual: "IIII" },
        FqRecord { name: "x2", seq: "TGCA", qual: "IIII" },
    ];
    let input = write_fastq(dir.path(), "input.fq", &reads);
    let output = dir.path().join("output.fq.gz");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "3",
        "--threads",
        "2",
    ]);

    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("read name mismatch"), "unexpected error: {stderr}");
}

#[test]
fn test_fastq_names_with_suffixes() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());

    let kraken = write_kraken_output(dir.path(), &["C\tr1\t3\t4\t3:1", "C\tr2\t3\t4\t3:1"]);
    let reads_r1 = vec![
        FqRecord { name: "r1/1", seq: "ACGT", qual: "IIII" },
        FqRecord { name: "r2/1", seq: "TGCA", qual: "IIII" },
    ];
    let reads_r2 = vec![
        FqRecord { name: "r1/2", seq: "TTTT", qual: "IIII" },
        FqRecord { name: "r2/2", seq: "AAAA", qual: "IIII" },
    ];
    let in1 = write_fastq(dir.path(), "r1.fq", &reads_r1);
    let in2 = write_fastq(dir.path(), "r2.fq", &reads_r2);
    let out1 = dir.path().join("out1.fq.gz");
    let out2 = dir.path().join("out2.fq.gz");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        in1.to_str().unwrap(),
        in2.to_str().unwrap(),
        "-o",
        out1.to_str().unwrap(),
        out2.to_str().unwrap(),
        "-t",
        "3",
        "--threads",
        "2",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));

    let recs1 = read_output_fastq(&out1);
    let recs2 = read_output_fastq(&out2);
    assert_eq!(recs1.len(), 2);
    assert_eq!(recs2.len(), 2);
    // Headers should preserve the original /1 /2 suffixes
    assert_eq!(recs1[0].0, "r1/1");
    assert_eq!(recs2[0].0, "r1/2");
}

#[test]
fn test_empty_inputs_single_end() {
    let dir = TempDir::new().unwrap();
    let report = write_empty_report(dir.path());
    let input = write_fastq(dir.path(), "input.fq", &[]);
    let output = dir.path().join("output.fq.gz");
    // kraken output file intentionally does not exist
    let kraken = dir.path().join("nonexistent_kraken.txt");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "1",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    assert!(output.exists());
    let records = read_output_fastq(&output);
    assert_eq!(records.len(), 0);
}

#[test]
fn test_empty_inputs_paired_end() {
    let dir = TempDir::new().unwrap();
    let report = write_empty_report(dir.path());
    let in1 = write_fastq(dir.path(), "r1.fq", &[]);
    let in2 = write_fastq(dir.path(), "r2.fq", &[]);
    let out1 = dir.path().join("out1.fq.gz");
    let out2 = dir.path().join("out2.fq.gz");
    let kraken = dir.path().join("nonexistent_kraken.txt");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        in1.to_str().unwrap(),
        in2.to_str().unwrap(),
        "-o",
        out1.to_str().unwrap(),
        out2.to_str().unwrap(),
        "-t",
        "1",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    assert_eq!(read_output_fastq(&out1).len(), 0);
    assert_eq!(read_output_fastq(&out2).len(), 0);
}

#[test]
fn test_empty_report_nonempty_fastq() {
    let dir = TempDir::new().unwrap();
    let report = write_empty_report(dir.path());
    let input = write_fastq(dir.path(), "input.fq", &make_reads());
    let output = dir.path().join("output.fq.gz");
    let kraken = dir.path().join("nonexistent_kraken.txt");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "1",
    ]);

    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("report is empty but FASTQ input") && stderr.contains("inconsistent"),
        "unexpected stderr: {stderr}"
    );
}

#[test]
fn test_nonempty_report_missing_kraken_output() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let input = write_fastq(dir.path(), "input.fq", &make_reads());
    let output = dir.path().join("output.fq.gz");
    let kraken = dir.path().join("nonexistent_kraken.txt");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "3",
    ]);

    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("failed to open kraken output"), "unexpected stderr: {stderr}");
}

#[test]
fn test_non_gz_extension_writes_uncompressed_output() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let input = write_fastq(dir.path(), "input.fq", &make_reads());
    let output = dir.path().join("output.fq");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "3",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    assert!(!is_gzip_file(&output));
    let records = parse_plain_fastq(&std::fs::read(&output).unwrap());
    assert_eq!(
        records,
        vec![("r1".to_string(), "ACGT".to_string()), ("r5".to_string(), "GGGG".to_string())]
    );
}

#[test]
fn test_bgz_extension_writes_bgzf_output() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let input = write_fastq(dir.path(), "input.fq", &make_reads());
    let output = dir.path().join("output.fq.bgz");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "3",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    assert!(is_gzip_file(&output));
    assert_eq!(read_output_fastq(&output).len(), 2);
}

#[test]
fn test_single_end_output_to_stdout() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let input = write_fastq(dir.path(), "input.fq", &make_reads());

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        "-",
        "-t",
        "3",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    let records = parse_plain_fastq(&result.stdout);
    assert_eq!(
        records,
        vec![("r1".to_string(), "ACGT".to_string()), ("r5".to_string(), "GGGG".to_string())]
    );
}

#[test]
fn test_paired_end_single_output_interleaves_to_stdout() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let in1 = write_fastq(dir.path(), "r1.fq", &make_reads());
    let in2 = write_fastq(dir.path(), "r2.fq", &make_reads_r2());

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        in1.to_str().unwrap(),
        in2.to_str().unwrap(),
        "-o",
        "-",
        "-t",
        "3",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    let records = parse_plain_fastq(&result.stdout);
    assert_eq!(
        records,
        vec![
            ("r1".to_string(), "ACGT".to_string()),
            ("r1".to_string(), "TTTT".to_string()),
            ("r5".to_string(), "GGGG".to_string()),
            ("r5".to_string(), "TTTT".to_string()),
        ]
    );
}

#[test]
fn test_paired_end_single_output_interleaves_to_bgzf_file() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let in1 = write_fastq(dir.path(), "r1.fq", &make_reads());
    let in2 = write_fastq(dir.path(), "r2.fq", &make_reads_r2());
    let output = dir.path().join("interleaved.fq.gz");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        in1.to_str().unwrap(),
        in2.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-t",
        "5",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    let records = read_output_fastq(&output);
    assert_eq!(
        records,
        vec![("r3".to_string(), "AAAA".to_string()), ("r3".to_string(), "CCCC".to_string())]
    );
}

#[test]
fn test_paired_end_mixed_uncompressed_and_bgzf_outputs() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let in1 = write_fastq(dir.path(), "r1.fq", &make_reads());
    let in2 = write_fastq(dir.path(), "r2.fq", &make_reads_r2());
    let out1 = dir.path().join("out1.fq");
    let out2 = dir.path().join("out2.fq.gz");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        in1.to_str().unwrap(),
        in2.to_str().unwrap(),
        "-o",
        out1.to_str().unwrap(),
        out2.to_str().unwrap(),
        "-t",
        "5",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    assert!(!is_gzip_file(&out1));
    assert!(is_gzip_file(&out2));
    assert_eq!(
        parse_plain_fastq(&std::fs::read(&out1).unwrap()),
        vec![("r3".to_string(), "AAAA".to_string())]
    );
    assert_eq!(read_output_fastq(&out2), vec![("r3".to_string(), "CCCC".to_string())]);
}

#[test]
fn test_stdout_for_both_outputs_errors() {
    let dir = TempDir::new().unwrap();
    let report = write_report(dir.path(), &standard_report_lines());
    let kraken = write_kraken_output(dir.path(), &make_kraken_lines());
    let in1 = write_fastq(dir.path(), "r1.fq", &make_reads());
    let in2 = write_fastq(dir.path(), "r2.fq", &make_reads_r2());

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        in1.to_str().unwrap(),
        in2.to_str().unwrap(),
        "-o",
        "-",
        "-",
        "-t",
        "3",
    ]);

    assert!(!result.status.success());
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("only one output may be written to stdout"), "unexpected: {stderr}");
}

#[test]
fn test_empty_inputs_to_stdout_writes_nothing() {
    let dir = TempDir::new().unwrap();
    let report = write_empty_report(dir.path());
    let input = write_fastq(dir.path(), "input.fq", &[]);
    let kraken = dir.path().join("nonexistent_kraken.txt");

    let result = run_filter(&[
        "-r",
        report.to_str().unwrap(),
        "-k",
        kraken.to_str().unwrap(),
        "-i",
        input.to_str().unwrap(),
        "-o",
        "-",
        "-t",
        "1",
    ]);

    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));
    assert!(result.stdout.is_empty());
}

#[test]
fn test_stdout_closed_early_by_reader_exits_successfully() {
    // Enough reads that the output far exceeds the pipe and write buffers, so writes
    // are still in progress when the reader goes away
    let num_reads = 20_000;
    let dir = TempDir::new().unwrap();
    let report = write_report(
        dir.path(),
        &[&format!("100.00\t{num_reads}\t{num_reads}\tS\t3\tEscherichia coli")],
    );
    let kraken_path = dir.path().join("kraken_output.txt");
    let fastq_path = dir.path().join("input.fq");
    let mut kraken = File::create(&kraken_path).unwrap();
    let mut fastq = File::create(&fastq_path).unwrap();
    for i in 0..num_reads {
        writeln!(kraken, "C\tread{i}\t3\t100\t3:1").unwrap();
        writeln!(fastq, "@read{i}\n{}\n+\n{}", "A".repeat(100), "I".repeat(100)).unwrap();
    }
    drop((kraken, fastq));

    let mut child = Command::new(k2tools_bin())
        .arg("filter")
        .args(["-r", report.to_str().unwrap()])
        .args(["-k", kraken_path.to_str().unwrap()])
        .args(["-i", fastq_path.to_str().unwrap()])
        .args(["-o", "-", "-t", "3"])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut stdout = child.stdout.take().unwrap();
    let mut first_bytes = [0_u8; 1024];
    stdout.read_exact(&mut first_bytes).unwrap();
    drop(stdout);
    let result = child.wait_with_output().unwrap();

    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(result.status.success(), "stderr: {stderr}");
    assert!(!stderr.contains("ERROR"), "unexpected error output: {stderr}");
    assert!(first_bytes.starts_with(b"@read0\n"));
}
