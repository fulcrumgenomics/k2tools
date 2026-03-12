pub mod kraken_report;
/// Parsing and representation of kraken2 report files.
///
/// Kraken2 reports are tab-delimited files with one row per taxon, encoding a taxonomy
/// tree via DFS order and indentation. This module supports both the standard 6-column
/// format and the extended 8-column format (with minimizer data).
pub mod rank;
pub mod row;

pub use kraken_report::KrakenReport;
pub use rank::{Rank, TaxonomicRank};
pub use row::ReportRow;
