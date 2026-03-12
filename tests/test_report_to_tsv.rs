use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::TempDir;

/// Returns the path to the k2tools binary.
fn k2tools_bin() -> PathBuf {
    let mut path = std::env::current_exe().unwrap();
    path.pop(); // remove test binary name
    path.pop(); // remove 'deps'
    path.push("k2tools");
    path
}

/// Runs `k2tools report-to-tsv` with the given arguments and returns the output.
fn run_report_to_tsv(args: &[&str]) -> std::process::Output {
    Command::new(k2tools_bin()).arg("report-to-tsv").args(args).output().unwrap()
}

/// Builds a standard (6-column) report line.
fn standard_line(
    pct: f64,
    clade: u64,
    direct: u64,
    rank: &str,
    taxid: u64,
    name: &str,
    depth: usize,
) -> String {
    let indent = " ".repeat(depth * 2);
    format!("{pct:.2}\t{clade}\t{direct}\t{rank}\t{taxid}\t{indent}{name}")
}

/// Builds an extended (8-column) report line with minimizer columns.
#[allow(clippy::too_many_arguments)]
fn extended_line(
    pct: f64,
    clade: u64,
    direct: u64,
    minimizers: u64,
    distinct: u64,
    rank: &str,
    taxid: u64,
    name: &str,
    depth: usize,
) -> String {
    let indent = " ".repeat(depth * 2);
    format!(
        "{pct:.2}\t{clade}\t{direct}\t{minimizers}\t{distinct}\t{rank}\t{taxid}\t{indent}{name}"
    )
}

/// Writes a report file and returns its path.
fn write_report(dir: &Path, lines: &[String]) -> PathBuf {
    let path = dir.join("report.txt");
    let content = lines.join("\n") + "\n";
    std::fs::write(&path, content).unwrap();
    path
}

/// Parses TSV output into a header and rows of field vectors.
fn parse_tsv(output: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let mut lines = output.lines();
    let header: Vec<String> =
        lines.next().expect("expected header line").split('\t').map(String::from).collect();

    let rows: Vec<Vec<String>> =
        lines.map(|line| line.split('\t').map(String::from).collect()).collect();

    (header, rows)
}

#[test]
fn test_standard_report_to_tsv() {
    let dir = TempDir::new().unwrap();
    let report_path = write_report(
        dir.path(),
        &[
            standard_line(10.0, 100, 100, "U", 0, "unclassified", 0),
            standard_line(90.0, 900, 5, "R", 1, "root", 0),
            standard_line(60.0, 600, 10, "D", 2, "Bacteria", 1),
            standard_line(50.0, 500, 500, "S", 3, "Escherichia coli", 2),
            standard_line(30.0, 300, 10, "D", 4, "Eukaryota", 1),
            standard_line(20.0, 200, 200, "S", 5, "Homo sapiens", 2),
        ],
    );

    let output = run_report_to_tsv(&["-r", report_path.to_str().unwrap()]);
    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));

    let stdout = String::from_utf8(output.stdout).unwrap();
    let (header, rows) = parse_tsv(&stdout);

    // Verify expected header columns
    assert_eq!(header.len(), 14);
    assert_eq!(header[0], "tax_id");
    assert_eq!(header[1], "name");
    assert_eq!(header[2], "rank");
    assert_eq!(header[3], "level");
    assert_eq!(header[4], "parent_tax_id");
    assert_eq!(header[5], "parent_rank");
    assert_eq!(header[6], "clade_count");
    assert_eq!(header[7], "direct_count");
    assert_eq!(header[8], "descendant_count");
    assert_eq!(header[9], "frac_clade");
    assert_eq!(header[10], "frac_direct");
    assert_eq!(header[11], "frac_descendant");
    assert_eq!(header[12], "minimizer_count");
    assert_eq!(header[13], "distinct_minimizer_count");

    // Verify row count
    assert_eq!(rows.len(), 6);

    // All rows have same column count as header
    for row in &rows {
        assert_eq!(row.len(), header.len());
    }

    // Check unclassified row
    assert_eq!(rows[0][0], "0");
    assert_eq!(rows[0][1], "unclassified");
    assert_eq!(rows[0][2], "U");
    assert_eq!(rows[0][3], "0"); // level
    assert_eq!(rows[0][4], ""); // no parent
    assert_eq!(rows[0][5], ""); // no parent
    assert_eq!(rows[0][6], "100"); // clade_count
    assert_eq!(rows[0][7], "100"); // direct_count
    assert_eq!(rows[0][8], "0"); // descendant_count

    // Check E. coli row
    assert_eq!(rows[3][0], "3");
    assert_eq!(rows[3][1], "Escherichia coli");
    assert_eq!(rows[3][2], "S");
    assert_eq!(rows[3][3], "2"); // level (child of Bacteria, grandchild of root)
    assert_eq!(rows[3][4], "2"); // parent is Bacteria
    assert_eq!(rows[3][5], "D"); // parent rank
    assert_eq!(rows[3][6], "500"); // clade_count
    assert_eq!(rows[3][7], "500"); // direct_count
    assert_eq!(rows[3][8], "0"); // descendant_count

    // Minimizer columns should be empty for standard report
    assert_eq!(rows[0][12], "");
    assert_eq!(rows[0][13], "");

    // Check a fraction value: E. coli clade=500/1000=0.5
    let frac: f64 = rows[3][9].parse().unwrap();
    assert!((frac - 0.5).abs() < 1e-9);
}

#[test]
fn test_output_to_file() {
    let dir = TempDir::new().unwrap();
    let report_path = write_report(
        dir.path(),
        &[
            standard_line(50.0, 5, 5, "U", 0, "unclassified", 0),
            standard_line(50.0, 5, 5, "R", 1, "root", 0),
        ],
    );
    let output_path = dir.path().join("output.tsv");

    let result = run_report_to_tsv(&[
        "-r",
        report_path.to_str().unwrap(),
        "-o",
        output_path.to_str().unwrap(),
    ]);
    assert!(result.status.success(), "stderr: {}", String::from_utf8_lossy(&result.stderr));

    let content = std::fs::read_to_string(&output_path).unwrap();
    let (header, rows) = parse_tsv(&content);
    assert_eq!(header.len(), 14);
    assert_eq!(rows.len(), 2);
}

#[test]
fn test_empty_report_produces_header_only() {
    let dir = TempDir::new().unwrap();
    let report_path = dir.path().join("report.txt");
    std::fs::write(&report_path, "").unwrap();

    let output = run_report_to_tsv(&["-r", report_path.to_str().unwrap()]);
    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));

    let stdout = String::from_utf8(output.stdout).unwrap();
    let (header, rows) = parse_tsv(&stdout);
    assert_eq!(header.len(), 14);
    assert!(rows.is_empty());
}

#[test]
fn test_extended_report_with_minimizers() {
    let dir = TempDir::new().unwrap();
    let report_path = write_report(
        dir.path(),
        &[
            extended_line(10.0, 100, 100, 0, 0, "U", 0, "unclassified", 0),
            extended_line(90.0, 900, 5, 500, 400, "R", 1, "root", 0),
            extended_line(60.0, 600, 10, 300, 250, "D", 2, "Bacteria", 1),
            extended_line(50.0, 500, 500, 200, 150, "S", 3, "Escherichia coli", 2),
        ],
    );

    let output = run_report_to_tsv(&["-r", report_path.to_str().unwrap()]);
    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));

    let stdout = String::from_utf8(output.stdout).unwrap();
    let (_, rows) = parse_tsv(&stdout);

    assert_eq!(rows.len(), 4);

    // Bacteria: minimizers=300, distinct=250
    assert_eq!(rows[2][12], "300");
    assert_eq!(rows[2][13], "250");

    // Unclassified: minimizers=0, distinct=0
    assert_eq!(rows[0][12], "0");
    assert_eq!(rows[0][13], "0");
}
