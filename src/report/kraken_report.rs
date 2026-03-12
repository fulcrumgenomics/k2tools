use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;
use std::str::FromStr;

use anyhow::{Context, Result};

use super::rank::TaxonomicRank;
use super::row::ReportRow;

/// The number of tab-delimited columns in a standard kraken2 report.
const STANDARD_COLUMNS: usize = 6;

/// The number of tab-delimited columns in an extended kraken2 report (with minimizer data).
const EXTENDED_COLUMNS: usize = 8;

/// A parsed kraken2 report, providing access to rows in DFS order, parent-child
/// relationships, and O(1) taxon ID lookups.
///
/// Reports are immutable after construction. Create one via [`KrakenReport::from_path`]
/// or [`KrakenReport::from_reader`].
pub struct KrakenReport {
    rows: Vec<ReportRow>,
    parents: Vec<Option<usize>>,
    taxon_id_to_index: HashMap<u64, usize>,
    has_minimizer_data: bool,
}

impl KrakenReport {
    /// Parses a kraken2 report from a file path.
    ///
    /// # Errors
    /// Returns an error if the file cannot be opened or the report content is malformed.
    pub fn from_path(path: &Path) -> Result<Self> {
        let file = std::fs::File::open(path)
            .with_context(|| format!("failed to open report file: {}", path.display()))?;
        let reader = std::io::BufReader::new(file);
        Self::from_reader(reader)
            .with_context(|| format!("failed to parse report file: {}", path.display()))
    }

    /// Parses a kraken2 report from any buffered reader.
    ///
    /// Auto-detects standard (6-column) vs extended (8-column) format from the first
    /// non-empty line.
    ///
    /// # Errors
    /// Returns an error if the input is empty, has an unsupported column count, contains
    /// inconsistent column counts, has non-numeric fields, duplicate taxon IDs, or
    /// odd leading-space counts in names.
    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self> {
        let lines: Vec<String> = reader
            .lines()
            .collect::<std::io::Result<Vec<_>>>()
            .context("failed to read report lines")?;

        let non_empty: Vec<&str> =
            lines.iter().map(String::as_str).filter(|l| !l.is_empty()).collect();
        if non_empty.is_empty() {
            return Ok(Self {
                rows: Vec::new(),
                parents: Vec::new(),
                taxon_id_to_index: HashMap::new(),
                has_minimizer_data: false,
            });
        }

        // Detect format from first line
        let expected_columns = detect_column_count(non_empty[0])?;
        let has_minimizer_data = expected_columns == EXTENDED_COLUMNS;

        let mut rows = Vec::with_capacity(non_empty.len());
        for (line_num, line) in non_empty.iter().enumerate() {
            let row = parse_row(line, expected_columns, has_minimizer_data)
                .with_context(|| format!("line {}", line_num + 1))?;
            rows.push(row);
        }

        let parents = build_parents(&rows);
        let taxon_id_to_index = build_taxon_index(&rows)?;

        Ok(Self { rows, parents, taxon_id_to_index, has_minimizer_data })
    }

    /// Returns all rows in DFS order.
    #[must_use]
    pub fn rows(&self) -> &[ReportRow] {
        &self.rows
    }

    /// Returns the number of rows in the report.
    #[must_use]
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    /// Returns `true` if the report contains no rows.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// Returns the row at the given index.
    ///
    /// # Panics
    /// Panics if `index >= self.len()`.
    #[must_use]
    pub fn row(&self, index: usize) -> &ReportRow {
        &self.rows[index]
    }

    /// Looks up a row by taxon ID, returning `None` if the ID is not present.
    #[must_use]
    pub fn get_by_taxon_id(&self, taxon_id: u64) -> Option<&ReportRow> {
        self.taxon_id_to_index.get(&taxon_id).map(|&i| &self.rows[i])
    }

    /// Returns the row index for a taxon ID, or `None` if not present.
    #[must_use]
    pub fn index_of_taxon_id(&self, taxon_id: u64) -> Option<usize> {
        self.taxon_id_to_index.get(&taxon_id).copied()
    }

    /// Returns the parent row, or `None` if the row is at the root level.
    ///
    /// # Panics
    /// Panics if `index >= self.len()`.
    #[must_use]
    pub fn parent(&self, index: usize) -> Option<&ReportRow> {
        self.parents[index].map(|pi| &self.rows[pi])
    }

    /// Returns the parent's row index, or `None` if the row is at the root level.
    ///
    /// # Panics
    /// Panics if `index >= self.len()`.
    #[must_use]
    pub fn parent_index(&self, index: usize) -> Option<usize> {
        self.parents[index]
    }

    /// Returns the indices of all descendants of the row at `index`.
    ///
    /// Leverages the DFS ordering: all descendants of row `i` are the contiguous
    /// run of rows after `i` with depth greater than `rows[i].depth()`.
    ///
    /// # Panics
    /// Panics if `index >= self.len()`.
    #[must_use]
    pub fn descendants(&self, index: usize) -> Vec<usize> {
        let parent_depth = self.rows[index].depth();
        let mut result = Vec::new();
        for i in (index + 1)..self.rows.len() {
            if self.rows[i].depth() <= parent_depth {
                break;
            }
            result.push(i);
        }
        result
    }

    /// Returns the indices of the direct children of the row at `index`.
    ///
    /// Computed on-demand by scanning forward in DFS order.
    ///
    /// # Panics
    /// Panics if `index >= self.len()`.
    #[must_use]
    pub fn children(&self, index: usize) -> Vec<usize> {
        let parent_depth = self.rows[index].depth();
        let mut result = Vec::new();
        for i in (index + 1)..self.rows.len() {
            let d = self.rows[i].depth();
            if d <= parent_depth {
                break;
            }
            if d == parent_depth + 1 {
                result.push(i);
            }
        }
        result
    }

    /// Returns the total number of sequences (classified + unclassified).
    ///
    /// Computed as the sum of clade counts for all depth-0 rows.
    #[must_use]
    pub fn total_sequences(&self) -> u64 {
        self.rows.iter().filter(|r| r.depth() == 0).map(ReportRow::clade_count).sum()
    }

    /// Returns `true` if this report contains minimizer data (extended 8-column format).
    #[must_use]
    pub fn has_minimizer_data(&self) -> bool {
        self.has_minimizer_data
    }
}

/// Detects the column count from a tab-delimited line.
fn detect_column_count(line: &str) -> Result<usize> {
    let count = line.split('\t').count();
    anyhow::ensure!(
        count == STANDARD_COLUMNS || count == EXTENDED_COLUMNS,
        "expected {STANDARD_COLUMNS} or {EXTENDED_COLUMNS} columns, found {count}"
    );
    Ok(count)
}

/// Parses a single tab-delimited report line into a `ReportRow`.
fn parse_row(line: &str, expected_columns: usize, has_minimizer_data: bool) -> Result<ReportRow> {
    let fields: Vec<&str> = line.split('\t').collect();
    anyhow::ensure!(
        fields.len() == expected_columns,
        "expected {expected_columns} columns, found {}",
        fields.len()
    );

    let percentage: f64 =
        fields[0].trim().parse().with_context(|| format!("invalid percentage: {:?}", fields[0]))?;
    let clade_count: u64 = fields[1]
        .trim()
        .parse()
        .with_context(|| format!("invalid clade count: {:?}", fields[1]))?;
    let direct_count: u64 = fields[2]
        .trim()
        .parse()
        .with_context(|| format!("invalid direct count: {:?}", fields[2]))?;

    // Column indices shift based on format
    let (minimizer_count, distinct_minimizer_count, rank_col, taxon_col, name_col) =
        if has_minimizer_data {
            let mc: u64 = fields[3]
                .trim()
                .parse()
                .with_context(|| format!("invalid minimizer count: {:?}", fields[3]))?;
            let dmc: u64 = fields[4]
                .trim()
                .parse()
                .with_context(|| format!("invalid distinct minimizer count: {:?}", fields[4]))?;
            (Some(mc), Some(dmc), 5, 6, 7)
        } else {
            (None, None, 3, 4, 5)
        };

    let taxonomic_rank = TaxonomicRank::from_str(fields[rank_col].trim())
        .with_context(|| format!("invalid rank code: {:?}", fields[rank_col]))?;
    let taxon_id: u64 = fields[taxon_col]
        .trim()
        .parse()
        .with_context(|| format!("invalid taxon ID: {:?}", fields[taxon_col]))?;

    // Derive depth from leading spaces in the name column
    let name_raw = fields[name_col];
    let leading_spaces = name_raw.len() - name_raw.trim_start_matches(' ').len();
    anyhow::ensure!(
        leading_spaces % 2 == 0,
        "name has odd number of leading spaces ({leading_spaces}): {name_raw:?}"
    );
    let depth = leading_spaces / 2;
    let name = name_raw.trim().to_string();

    Ok(ReportRow::new(
        percentage,
        clade_count,
        direct_count,
        minimizer_count,
        distinct_minimizer_count,
        taxonomic_rank,
        taxon_id,
        name,
        depth,
    ))
}

/// Builds the parent index array using a stack-based DFS reconstruction.
fn build_parents(rows: &[ReportRow]) -> Vec<Option<usize>> {
    let mut parents = vec![None; rows.len()];
    // Stack of (row_index, depth)
    let mut stack: Vec<(usize, usize)> = Vec::new();

    for (i, row) in rows.iter().enumerate() {
        let d = row.depth();
        while let Some(&(_, top_depth)) = stack.last() {
            if top_depth >= d {
                stack.pop();
            } else {
                break;
            }
        }
        parents[i] = stack.last().map(|&(idx, _)| idx);
        stack.push((i, d));
    }

    parents
}

/// Builds a taxon ID to row index map, erroring on duplicates.
fn build_taxon_index(rows: &[ReportRow]) -> Result<HashMap<u64, usize>> {
    let mut map = HashMap::with_capacity(rows.len());
    for (i, row) in rows.iter().enumerate() {
        let prev = map.insert(row.taxon_id(), i);
        anyhow::ensure!(
            prev.is_none(),
            "duplicate taxon ID {} at rows {} and {}",
            row.taxon_id(),
            prev.unwrap_or(0),
            i
        );
    }
    Ok(map)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::rank::Rank;

    /// Builds a standard 6-column report line.
    fn standard_line(
        pct: f64,
        clade: u64,
        direct: u64,
        rank: &str,
        taxid: u64,
        name: &str,
        depth: usize,
    ) -> String {
        let indent = "  ".repeat(depth);
        format!("{pct:6.2}\t{clade}\t{direct}\t{rank}\t{taxid}\t{indent}{name}")
    }

    /// Builds an extended 8-column report line.
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
        let indent = "  ".repeat(depth);
        format!(
            "{pct:6.2}\t{clade}\t{direct}\t{minimizers}\t{distinct}\t{rank}\t{taxid}\t{indent}{name}"
        )
    }

    /// Builds a minimal standard report string.
    fn minimal_standard_report() -> String {
        [
            standard_line(10.0, 100, 100, "U", 0, "unclassified", 0),
            standard_line(90.0, 900, 5, "R", 1, "root", 0),
            standard_line(90.0, 895, 895, "S", 9606, "Homo sapiens", 1),
        ]
        .join("\n")
    }

    fn parse(report: &str) -> KrakenReport {
        KrakenReport::from_reader(report.as_bytes()).unwrap()
    }

    #[test]
    fn test_parse_minimal_standard_report() {
        let kr = parse(&minimal_standard_report());
        assert_eq!(kr.len(), 3);
        assert!(!kr.has_minimizer_data());

        assert_eq!(kr.row(0).name(), "unclassified");
        assert_eq!(kr.row(0).taxon_id(), 0);
        assert_eq!(kr.row(0).depth(), 0);
        assert_eq!(kr.row(0).clade_count(), 100);

        assert_eq!(kr.row(1).name(), "root");
        assert_eq!(kr.row(1).taxon_id(), 1);

        assert_eq!(kr.row(2).name(), "Homo sapiens");
        assert_eq!(kr.row(2).taxon_id(), 9606);
        assert_eq!(kr.row(2).depth(), 1);
    }

    #[test]
    fn test_parse_extended_report() {
        let report = [
            extended_line(10.0, 100, 100, 0, 0, "U", 0, "unclassified", 0),
            extended_line(90.0, 900, 5, 5000, 3000, "R", 1, "root", 0),
            extended_line(90.0, 895, 895, 4500, 2800, "S", 9606, "Homo sapiens", 1),
        ]
        .join("\n");

        let kr = parse(&report);
        assert!(kr.has_minimizer_data());
        assert_eq!(kr.row(1).minimizer_count(), Some(5000));
        assert_eq!(kr.row(1).distinct_minimizer_count(), Some(3000));
    }

    #[test]
    fn test_standard_report_no_minimizers() {
        let kr = parse(&minimal_standard_report());
        assert_eq!(kr.row(0).minimizer_count(), None);
        assert_eq!(kr.row(0).distinct_minimizer_count(), None);
    }

    #[test]
    fn test_parent_child_multi_level() {
        let report = [
            standard_line(100.0, 1000, 10, "R", 1, "root", 0),
            standard_line(80.0, 800, 5, "D", 2, "Bacteria", 1),
            standard_line(60.0, 600, 3, "P", 3, "Proteobacteria", 2),
            standard_line(40.0, 400, 400, "C", 4, "Gammaproteobacteria", 3),
        ]
        .join("\n");

        let kr = parse(&report);

        // root has no parent
        assert!(kr.parent(0).is_none());
        assert_eq!(kr.parent_index(0), None);

        // Bacteria -> root
        assert_eq!(kr.parent_index(1), Some(0));
        assert_eq!(kr.parent(1).unwrap().name(), "root");

        // Proteobacteria -> Bacteria
        assert_eq!(kr.parent_index(2), Some(1));

        // Gammaproteobacteria -> Proteobacteria
        assert_eq!(kr.parent_index(3), Some(2));
    }

    #[test]
    fn test_children_multiple() {
        let report = [
            standard_line(100.0, 1000, 10, "R", 1, "root", 0),
            standard_line(60.0, 600, 5, "D", 2, "Bacteria", 1),
            standard_line(40.0, 400, 400, "P", 3, "Proteobacteria", 2),
            standard_line(40.0, 400, 5, "D", 4, "Eukaryota", 1),
            standard_line(30.0, 300, 300, "P", 5, "Chordata", 2),
        ]
        .join("\n");

        let kr = parse(&report);
        let root_children = kr.children(0);
        assert_eq!(root_children, vec![1, 3]);

        let bacteria_children = kr.children(1);
        assert_eq!(bacteria_children, vec![2]);

        let eukaryota_children = kr.children(3);
        assert_eq!(eukaryota_children, vec![4]);

        // Leaf has no children
        assert!(kr.children(2).is_empty());
    }

    #[test]
    fn test_taxon_id_lookup_hit() {
        let kr = parse(&minimal_standard_report());
        let row = kr.get_by_taxon_id(9606).unwrap();
        assert_eq!(row.name(), "Homo sapiens");
    }

    #[test]
    fn test_taxon_id_lookup_miss() {
        let kr = parse(&minimal_standard_report());
        assert!(kr.get_by_taxon_id(99999).is_none());
    }

    #[test]
    fn test_index_of_taxon_id() {
        let kr = parse(&minimal_standard_report());
        assert_eq!(kr.index_of_taxon_id(1), Some(1));
        assert_eq!(kr.index_of_taxon_id(99999), None);
    }

    #[test]
    fn test_total_sequences() {
        let kr = parse(&minimal_standard_report());
        // unclassified(100) + root(900)
        assert_eq!(kr.total_sequences(), 1000);
    }

    #[test]
    fn test_total_sequences_no_unclassified() {
        let report = [
            standard_line(100.0, 1000, 10, "R", 1, "root", 0),
            standard_line(100.0, 990, 990, "S", 9606, "Homo sapiens", 1),
        ]
        .join("\n");

        let kr = parse(&report);
        assert_eq!(kr.total_sequences(), 1000);
    }

    #[test]
    fn test_depth_derivation() {
        let report = [
            standard_line(100.0, 1000, 10, "R", 1, "root", 0),
            standard_line(80.0, 800, 5, "D", 2, "Bacteria", 1),
            standard_line(60.0, 600, 3, "P", 3, "Proteobacteria", 2),
            standard_line(40.0, 400, 400, "C", 4, "Gamma", 3),
        ]
        .join("\n");

        let kr = parse(&report);
        assert_eq!(kr.row(0).depth(), 0);
        assert_eq!(kr.row(1).depth(), 1);
        assert_eq!(kr.row(2).depth(), 2);
        assert_eq!(kr.row(3).depth(), 3);
    }

    #[test]
    fn test_names_are_trimmed() {
        let kr = parse(&minimal_standard_report());
        assert_eq!(kr.row(2).name(), "Homo sapiens");
        assert!(!kr.row(2).name().starts_with(' '));
    }

    #[test]
    fn test_empty_input_returns_empty_report() {
        let kr = KrakenReport::from_reader("".as_bytes()).unwrap();
        assert!(kr.is_empty());
        assert_eq!(kr.len(), 0);
        assert_eq!(kr.total_sequences(), 0);
    }

    #[test]
    fn test_only_blank_lines_returns_empty_report() {
        let kr = KrakenReport::from_reader("\n\n\n".as_bytes()).unwrap();
        assert!(kr.is_empty());
        assert_eq!(kr.len(), 0);
    }

    #[test]
    fn test_error_inconsistent_columns() {
        let report = [
            standard_line(100.0, 1000, 10, "R", 1, "root", 0),
            "50.00\t500\t500\tS\t9606".to_string(),
        ]
        .join("\n");

        let result = KrakenReport::from_reader(report.as_bytes());
        assert!(result.is_err());
    }

    #[test]
    fn test_error_non_numeric_field() {
        let report = " 50.00\tabc\t100\tR\t1\troot";
        let result = KrakenReport::from_reader(report.as_bytes());
        assert!(result.is_err());
    }

    #[test]
    fn test_error_duplicate_taxon_ids() {
        let report = [
            standard_line(50.0, 500, 500, "R", 1, "root", 0),
            standard_line(50.0, 500, 500, "S", 1, "also root", 1),
        ]
        .join("\n");

        let result = KrakenReport::from_reader(report.as_bytes());
        assert!(result.is_err());
    }

    #[test]
    fn test_error_odd_leading_spaces() {
        let report = " 50.00\t500\t500\tS\t9606\t   Homo sapiens";
        let result = KrakenReport::from_reader(report.as_bytes());
        assert!(result.is_err());
    }

    #[test]
    fn test_error_unsupported_column_count() {
        let report = "a\tb\tc\td";
        let result = KrakenReport::from_reader(report.as_bytes());
        assert!(result.is_err());
    }

    #[test]
    fn test_non_standard_ranks_in_tree() {
        let report = [
            standard_line(100.0, 1000, 10, "R", 1, "root", 0),
            standard_line(80.0, 800, 5, "G", 2, "Staphylococcus", 1),
            standard_line(60.0, 600, 3, "G1", 3, "subgenus A", 2),
            standard_line(40.0, 400, 400, "S", 4, "Staphylococcus aureus", 3),
            standard_line(20.0, 200, 200, "S1", 5, "Staphylococcus aureus subsp. aureus", 4),
        ]
        .join("\n");

        let kr = parse(&report);
        assert_eq!(kr.len(), 5);

        let g1 = kr.row(2);
        assert_eq!(g1.taxonomic_rank().rank(), Rank::Genus);
        assert_eq!(g1.taxonomic_rank().depth(), Some(1));
        assert!(!g1.taxonomic_rank().is_standard());

        let s1 = kr.row(4);
        assert_eq!(s1.taxonomic_rank().rank(), Rank::Species);
        assert_eq!(s1.taxonomic_rank().depth(), Some(1));

        // Parent chain: S1 -> S -> G1 -> G -> R
        assert_eq!(kr.parent_index(4), Some(3));
        assert_eq!(kr.parent_index(3), Some(2));
        assert_eq!(kr.parent_index(2), Some(1));
        assert_eq!(kr.parent_index(1), Some(0));
    }

    #[test]
    fn test_descendants_of_leaf() {
        let kr = parse(&minimal_standard_report());
        // Homo sapiens (index 2) is a leaf
        assert!(kr.descendants(2).is_empty());
    }

    #[test]
    fn test_descendants_of_root() {
        let report = [
            standard_line(100.0, 1000, 10, "R", 1, "root", 0),
            standard_line(80.0, 800, 5, "D", 2, "Bacteria", 1),
            standard_line(60.0, 600, 600, "P", 3, "Proteobacteria", 2),
        ]
        .join("\n");

        let kr = parse(&report);
        assert_eq!(kr.descendants(0), vec![1, 2]);
    }

    #[test]
    fn test_descendants_mid_level() {
        let report = [
            standard_line(100.0, 1000, 10, "R", 1, "root", 0),
            standard_line(60.0, 600, 5, "D", 2, "Bacteria", 1),
            standard_line(40.0, 400, 400, "P", 3, "Proteobacteria", 2),
            standard_line(40.0, 400, 5, "D", 4, "Eukaryota", 1),
            standard_line(30.0, 300, 300, "P", 5, "Chordata", 2),
        ]
        .join("\n");

        let kr = parse(&report);
        // Bacteria descendants = only Proteobacteria, not Eukaryota/Chordata
        assert_eq!(kr.descendants(1), vec![2]);
        // Eukaryota descendants = Chordata
        assert_eq!(kr.descendants(3), vec![4]);
    }

    #[test]
    fn test_is_empty() {
        let kr = parse(&minimal_standard_report());
        assert!(!kr.is_empty());
    }

    #[test]
    fn test_percentage_parsing() {
        let kr = parse(&minimal_standard_report());
        assert!((kr.row(0).percentage() - 10.0).abs() < f64::EPSILON);
        assert!((kr.row(1).percentage() - 90.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_blank_lines_are_skipped() {
        let report = format!(
            "\n{}\n\n{}\n{}\n",
            standard_line(10.0, 100, 100, "U", 0, "unclassified", 0),
            standard_line(90.0, 900, 5, "R", 1, "root", 0),
            standard_line(90.0, 895, 895, "S", 9606, "Homo sapiens", 1),
        );

        let kr = parse(&report);
        assert_eq!(kr.len(), 3);
    }
}
