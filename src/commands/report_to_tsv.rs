use std::fs::File;
use std::io::{BufWriter, Write, stdout};
use std::path::PathBuf;

use anyhow::Result;
use fgoxide::io::DelimFileWriter;
use serde::Serialize;

use crate::commands::command::Command;
use crate::report::KrakenReport;

/// Convert a kraken2 report file to a clean, header-bearing TSV.
///
/// Reads a kraken2 report (standard 6-column or extended 8-column with minimizer
/// data) and writes a tab-separated file with clearly named columns, derived
/// parent information, descendant counts, and fraction columns.
///
/// Output columns:
///
///   tax_id                   - NCBI taxonomy ID
///   name                     - Scientific name
///   rank                     - Taxonomic rank code (e.g. S, G, D1)
///   level                    - Depth in taxonomy (0 for root/unclassified)
///   parent_tax_id            - Parent taxon ID (empty for root/unclassified)
///   parent_rank              - Parent rank code (empty for root/unclassified)
///   clade_count              - Fragments in clade rooted at this taxon
///   direct_count             - Fragments assigned directly to this taxon
///   descendant_count         - clade_count minus direct_count
///   frac_clade               - clade_count / total_sequences
///   frac_direct              - direct_count / total_sequences
///   frac_descendant          - descendant_count / total_sequences
///   minimizer_count          - Minimizers in clade (empty if not in report)
///   distinct_minimizer_count - Distinct minimizers (empty if not in report)
///
/// Examples:
///
///   k2tools report-to-tsv -r kraken2_report.txt -o report.tsv
///   k2tools report-to-tsv -r kraken2_report.txt   # writes to stdout
#[derive(clap::Args)]
pub struct ReportToTsv {
    /// Path to the kraken2 report file.
    #[arg(short = 'r', long)]
    kraken_report: PathBuf,

    /// Output TSV file path. Defaults to stdout.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

/// One row of TSV output, assembled from the report before writing.
#[derive(Default, Serialize)]
struct TsvRow {
    tax_id: u64,
    name: String,
    rank: String,
    level: usize,
    parent_tax_id: String,
    parent_rank: String,
    clade_count: u64,
    direct_count: u64,
    descendant_count: u64,
    frac_clade: f64,
    frac_direct: f64,
    frac_descendant: f64,
    minimizer_count: String,
    distinct_minimizer_count: String,
}

/// Derives the TSV header line from `TsvRow`'s field names by serializing a default
/// row through a csv writer and extracting just the header.
fn tsv_header() -> String {
    let mut csv_writer =
        csv::WriterBuilder::new().delimiter(b'\t').has_headers(true).from_writer(Vec::new());
    csv_writer.serialize(TsvRow::default()).unwrap();
    csv_writer.flush().unwrap();
    let bytes = csv_writer.into_inner().unwrap();
    let text = String::from_utf8(bytes).unwrap();
    text.lines().next().unwrap().to_string()
}

impl Command for ReportToTsv {
    fn execute(&self) -> Result<()> {
        let report = KrakenReport::from_path(&self.kraken_report)?;
        let rows = build_tsv_rows(&report);

        let writer: BufWriter<Box<dyn Write + Send>> = match &self.output {
            Some(path) => {
                let file = File::create(path).map_err(|e| {
                    anyhow::anyhow!("failed to create output {}: {e}", path.display())
                })?;
                BufWriter::new(Box::new(file))
            }
            None => BufWriter::new(Box::new(stdout())),
        };

        if rows.is_empty() {
            // csv::Writer only emits headers on the first serialize call, so with
            // zero rows we derive the header from TsvRow's field names directly.
            let mut w = writer;
            writeln!(w, "{}", tsv_header())?;
            w.flush()?;
        } else {
            let mut tsv_writer = DelimFileWriter::new(writer, b'\t', true);
            tsv_writer.write_all(rows)?;
        }

        log::info!("Wrote {} rows to TSV.", report.len());
        Ok(())
    }
}

/// Builds `TsvRow` structs from every row in the report.
fn build_tsv_rows(report: &KrakenReport) -> Vec<TsvRow> {
    let total_sequences = report.total_sequences();
    let has_minimizer_data = report.has_minimizer_data();
    let mut tsv_rows = Vec::with_capacity(report.len());

    for (i, row) in report.rows().iter().enumerate() {
        let clade_count = row.clade_count();
        let direct_count = row.direct_count();
        let descendant_count = clade_count - direct_count;

        #[allow(clippy::cast_precision_loss)]
        let (frac_clade, frac_direct, frac_descendant) = if total_sequences > 0 {
            (
                clade_count as f64 / total_sequences as f64,
                direct_count as f64 / total_sequences as f64,
                descendant_count as f64 / total_sequences as f64,
            )
        } else {
            (0.0, 0.0, 0.0)
        };

        let (minimizer_count, distinct_minimizer_count) = if has_minimizer_data {
            let mc = row.minimizer_count().unwrap_or(0);
            let dmc = row.distinct_minimizer_count().unwrap_or(0);
            (mc.to_string(), dmc.to_string())
        } else {
            (String::new(), String::new())
        };

        let (parent_tax_id, parent_rank) = match report.parent(i) {
            Some(parent) => (parent.taxon_id().to_string(), parent.taxonomic_rank().to_string()),
            None => (String::new(), String::new()),
        };

        tsv_rows.push(TsvRow {
            tax_id: row.taxon_id(),
            name: row.name().to_string(),
            rank: row.taxonomic_rank().to_string(),
            level: row.depth(),
            parent_tax_id,
            parent_rank,
            clade_count,
            direct_count,
            descendant_count,
            frac_clade,
            frac_direct,
            frac_descendant,
            minimizer_count,
            distinct_minimizer_count,
        });
    }

    tsv_rows
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds a standard (6-column) report string from structured parameters.
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

    /// Builds an extended (8-column) report string with minimizer columns.
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

    fn parse(report: &str) -> KrakenReport {
        KrakenReport::from_reader(report.as_bytes()).unwrap()
    }

    fn make_standard_report() -> KrakenReport {
        parse(
            &[
                standard_line(10.0, 100, 100, "U", 0, "unclassified", 0),
                standard_line(90.0, 900, 5, "R", 1, "root", 0),
                standard_line(60.0, 600, 10, "D", 2, "Bacteria", 1),
                standard_line(50.0, 500, 500, "S", 3, "Escherichia coli", 2),
                standard_line(30.0, 300, 10, "D", 4, "Eukaryota", 1),
                standard_line(20.0, 200, 200, "S", 5, "Homo sapiens", 2),
            ]
            .join("\n"),
        )
    }

    fn make_extended_report() -> KrakenReport {
        parse(
            &[
                extended_line(10.0, 100, 100, 0, 0, "U", 0, "unclassified", 0),
                extended_line(90.0, 900, 5, 500, 400, "R", 1, "root", 0),
                extended_line(60.0, 600, 10, 300, 250, "D", 2, "Bacteria", 1),
                extended_line(50.0, 500, 500, 200, 150, "S", 3, "Escherichia coli", 2),
                extended_line(30.0, 300, 10, 200, 150, "D", 4, "Eukaryota", 1),
                extended_line(20.0, 200, 200, 100, 80, "S", 5, "Homo sapiens", 2),
            ]
            .join("\n"),
        )
    }

    /// Writes rows to a temp file using `DelimFile::write_tsv` and returns the TSV text.
    fn write_rows_to_string(rows: Vec<TsvRow>) -> String {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("out.tsv");
        let df = fgoxide::io::DelimFile::default();
        df.write_tsv(&path, rows).unwrap();
        std::fs::read_to_string(path).unwrap()
    }

    #[test]
    fn test_level_values() {
        let report = make_standard_report();
        let rows = build_tsv_rows(&report);

        assert_eq!(rows[0].level, 0); // unclassified
        assert_eq!(rows[1].level, 0); // root
        assert_eq!(rows[2].level, 1); // Bacteria (child of root)
        assert_eq!(rows[3].level, 2); // E. coli (child of Bacteria)
        assert_eq!(rows[4].level, 1); // Eukaryota (child of root)
        assert_eq!(rows[5].level, 2); // Homo sapiens (child of Eukaryota)
    }

    #[test]
    fn test_parent_fields_empty_for_root_and_unclassified() {
        let report = make_standard_report();
        let rows = build_tsv_rows(&report);

        // unclassified (depth 0, no parent)
        assert_eq!(rows[0].parent_tax_id, "");
        assert_eq!(rows[0].parent_rank, "");

        // root (depth 0, no parent)
        assert_eq!(rows[1].parent_tax_id, "");
        assert_eq!(rows[1].parent_rank, "");
    }

    #[test]
    fn test_parent_fields_populated_for_children() {
        let report = make_standard_report();
        let rows = build_tsv_rows(&report);

        // Bacteria -> parent is root (tax_id 1, rank R)
        assert_eq!(rows[2].parent_tax_id, "1");
        assert_eq!(rows[2].parent_rank, "R");

        // E. coli -> parent is Bacteria (tax_id 2, rank D)
        assert_eq!(rows[3].parent_tax_id, "2");
        assert_eq!(rows[3].parent_rank, "D");

        // Homo sapiens -> parent is Eukaryota (tax_id 4, rank D)
        assert_eq!(rows[5].parent_tax_id, "4");
        assert_eq!(rows[5].parent_rank, "D");
    }

    #[test]
    fn test_descendant_count() {
        let report = make_standard_report();
        let rows = build_tsv_rows(&report);

        // unclassified: clade=100, direct=100, desc=0
        assert_eq!(rows[0].descendant_count, 0);

        // root: clade=900, direct=5, desc=895
        assert_eq!(rows[1].descendant_count, 895);

        // Bacteria: clade=600, direct=10, desc=590
        assert_eq!(rows[2].descendant_count, 590);

        // E. coli: clade=500, direct=500, desc=0
        assert_eq!(rows[3].descendant_count, 0);
    }

    #[test]
    fn test_minimizer_columns_empty_for_standard_report() {
        let report = make_standard_report();
        let rows = build_tsv_rows(&report);

        for row in &rows {
            assert_eq!(row.minimizer_count, "");
            assert_eq!(row.distinct_minimizer_count, "");
        }
    }

    #[test]
    fn test_minimizer_columns_populated_for_extended_report() {
        let report = make_extended_report();
        let rows = build_tsv_rows(&report);

        // Root row: minimizers=500, distinct=400
        assert_eq!(rows[1].minimizer_count, "500");
        assert_eq!(rows[1].distinct_minimizer_count, "400");

        // Bacteria: minimizers=300, distinct=250
        assert_eq!(rows[2].minimizer_count, "300");
        assert_eq!(rows[2].distinct_minimizer_count, "250");

        // All minimizer count fields should be non-empty
        for row in &rows {
            assert!(!row.minimizer_count.is_empty());
            assert!(!row.distinct_minimizer_count.is_empty());
        }
    }

    #[test]
    fn test_fraction_calculations() {
        let report = make_standard_report();
        let rows = build_tsv_rows(&report);

        // total_sequences = 100 (unclassified) + 900 (root) = 1000
        // E. coli: clade=500, direct=500, descendant=0
        let ecoli = &rows[3];
        assert!((ecoli.frac_clade - 0.5).abs() < 1e-9);
        assert!((ecoli.frac_direct - 0.5).abs() < 1e-9);
        assert!((ecoli.frac_descendant - 0.0).abs() < 1e-9);

        // Bacteria: clade=600, direct=10, descendant=590
        let bacteria = &rows[2];
        assert!((bacteria.frac_clade - 0.6).abs() < 1e-9);
        assert!((bacteria.frac_direct - 0.01).abs() < 1e-9);
        assert!((bacteria.frac_descendant - 0.59).abs() < 1e-9);
    }

    #[test]
    fn test_fractions_zero_when_no_sequences() {
        let report = parse(
            &[
                standard_line(0.0, 0, 0, "U", 0, "unclassified", 0),
                standard_line(0.0, 0, 0, "R", 1, "root", 0),
            ]
            .join("\n"),
        );
        let rows = build_tsv_rows(&report);

        for row in &rows {
            assert!((row.frac_clade - 0.0).abs() < 1e-9);
            assert!((row.frac_direct - 0.0).abs() < 1e-9);
            assert!((row.frac_descendant - 0.0).abs() < 1e-9);
        }
    }

    #[test]
    fn test_basic_field_values() {
        let report = make_standard_report();
        let rows = build_tsv_rows(&report);

        assert_eq!(rows[0].tax_id, 0);
        assert_eq!(rows[0].name, "unclassified");
        assert_eq!(rows[0].rank, "U");
        assert_eq!(rows[0].clade_count, 100);
        assert_eq!(rows[0].direct_count, 100);

        assert_eq!(rows[3].tax_id, 3);
        assert_eq!(rows[3].name, "Escherichia coli");
        assert_eq!(rows[3].rank, "S");
    }

    #[test]
    fn test_empty_report_produces_no_rows() {
        let report = parse("");
        let rows = build_tsv_rows(&report);
        assert!(rows.is_empty());
    }

    #[test]
    fn test_write_tsv_header_and_rows() {
        let report = parse(
            &[
                standard_line(50.0, 5, 5, "U", 0, "unclassified", 0),
                standard_line(50.0, 5, 5, "R", 1, "root", 0),
            ]
            .join("\n"),
        );
        let rows = build_tsv_rows(&report);
        let text = write_rows_to_string(rows);
        let lines: Vec<&str> = text.lines().collect();

        // Header + 2 data rows
        assert_eq!(lines.len(), 3);

        // Header should match field names
        let header_cols: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(header_cols[0], "tax_id");
        assert_eq!(header_cols[3], "level");
        assert_eq!(header_cols[6], "clade_count");
        assert_eq!(header_cols[9], "frac_clade");
        assert_eq!(header_cols[12], "minimizer_count");
        assert_eq!(header_cols[13], "distinct_minimizer_count");
        assert_eq!(header_cols.len(), 14);

        // Verify column count consistency
        for line in &lines[1..] {
            assert_eq!(line.split('\t').count(), header_cols.len());
        }
    }
}
