use serde::{Deserialize, Serialize};

use super::rank::TaxonomicRank;

/// A single row from a kraken2 report file.
///
/// Each row represents one taxon and includes classification counts, the taxonomic rank,
/// and the taxon name. The `depth` field is derived from the leading whitespace indentation
/// in the original report and encodes the position in the taxonomy tree.
///
/// For extended reports (with minimizer data), `minimizer_count` and
/// `distinct_minimizer_count` are `Some`; for standard reports they are `None`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReportRow {
    percentage: f64,
    clade_count: u64,
    direct_count: u64,
    minimizer_count: Option<u64>,
    distinct_minimizer_count: Option<u64>,
    taxonomic_rank: TaxonomicRank,
    taxon_id: u64,
    name: String,
    depth: usize,
}

impl ReportRow {
    /// Creates a new `ReportRow`. Not public — rows are created only during report parsing.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        percentage: f64,
        clade_count: u64,
        direct_count: u64,
        minimizer_count: Option<u64>,
        distinct_minimizer_count: Option<u64>,
        taxonomic_rank: TaxonomicRank,
        taxon_id: u64,
        name: String,
        depth: usize,
    ) -> Self {
        Self {
            percentage,
            clade_count,
            direct_count,
            minimizer_count,
            distinct_minimizer_count,
            taxonomic_rank,
            taxon_id,
            name,
            depth,
        }
    }

    /// Percentage of fragments rooted at this taxon's clade.
    #[must_use]
    pub fn percentage(&self) -> f64 {
        self.percentage
    }

    /// Number of fragments in the clade rooted at this taxon.
    #[must_use]
    pub fn clade_count(&self) -> u64 {
        self.clade_count
    }

    /// Number of fragments assigned directly to this taxon.
    #[must_use]
    pub fn direct_count(&self) -> u64 {
        self.direct_count
    }

    /// Number of minimizers in the clade (extended reports only).
    #[must_use]
    pub fn minimizer_count(&self) -> Option<u64> {
        self.minimizer_count
    }

    /// Number of distinct minimizers in the clade (extended reports only).
    #[must_use]
    pub fn distinct_minimizer_count(&self) -> Option<u64> {
        self.distinct_minimizer_count
    }

    /// The taxonomic rank of this taxon.
    #[must_use]
    pub fn taxonomic_rank(&self) -> TaxonomicRank {
        self.taxonomic_rank
    }

    /// NCBI taxonomy ID (0 for unclassified).
    #[must_use]
    pub fn taxon_id(&self) -> u64 {
        self.taxon_id
    }

    /// Scientific name of the taxon (trimmed, no leading indentation).
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Depth in the taxonomy tree, derived from leading indentation (spaces / 2).
    #[must_use]
    pub fn depth(&self) -> usize {
        self.depth
    }
}
