use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use fgoxide::io::Io;
use fgoxide::iter::IntoChunkedReadAheadIterator;
use pooled_writer::bgzf::BgzfCompressor;
use pooled_writer::{Pool, PoolBuilder, PooledWriter};
use seq_io::fastq::{Error as FastqError, OwnedRecord, Reader as FastqReader, Record};

use crate::commands::command::Command;
use crate::kraken_output::{KrakenOutputReader, KrakenRecord};
use crate::progress::{ProgressLogger, format_count};
use crate::report::{KrakenReport, ReportRow};

/// Number of records per chunk sent through the read-ahead channel.
const READ_AHEAD_CHUNK_SIZE: usize = 1024;

/// Number of buffered chunks in the read-ahead channel.
const READ_AHEAD_NUM_CHUNKS: usize = 1024;

/// Buffer size used when opening input files for reading and uncompressed outputs for
/// writing.
const IO_BUFFER_SIZE: usize = 512 * 1024;

/// Output path that writes uncompressed FASTQ to stdout instead of a file.
const STDOUT_PATH: &str = "-";

/// Filter reads from FASTQ files based on kraken2 classification results
///
/// Extracts reads classified to one or more taxon IDs from FASTQ files, using the
/// kraken2 report (taxonomy tree) and per-read classification output. Supports both
/// single-end and paired-end reads, and writes bgzf-compressed or uncompressed FASTQ to
/// files or stdout.
///
/// # Required inputs
///
/// The command needs three pieces of data that must all come from the same kraken2 run:
///
/// - `--kraken-report` (`-r`): The kraken2 report file containing the taxonomy tree and
///   per-taxon read counts. This is used to resolve taxon IDs, expand descendants, and
///   estimate the expected number of matching reads.
/// - `--kraken-output` (`-k`): The per-read classification output from kraken2
///   (generated with `--output`). Each line maps a read name to a taxon ID.
/// - `--input` (`-i`): One FASTQ file for single-end data, or two for paired-end. Gzip
///   and bgzf compressed inputs are detected and handled automatically.
///
/// The kraken output and FASTQ file(s) must contain the same reads in the same order.
/// The command verifies read name agreement and will error if the files are mismatched
/// or have different numbers of records.
///
/// # Taxon selection
///
/// At least one of `--taxon-ids` or `--include-unclassified` must be specified.
///
/// - `--taxon-ids` (`-t`): One or more NCBI taxon IDs to extract. By default, only reads
///   classified directly to these exact taxon IDs are included.
/// - `--include-descendants` (`-d`): Expand each taxon ID to include all of its
///   descendants in the taxonomy tree. For example, specifying a genus-level taxon ID
///   with `-d` will also extract reads classified to any species or strain within that
///   genus.
/// - `--include-unclassified` (`-u`): Include reads that kraken2 could not classify
///   (taxon ID 0). Can be combined with `--taxon-ids` to extract both classified and
///   unclassified reads in a single pass.
/// - `--exclude-taxon-ids` (`-e`): Taxon IDs to remove from the selection built by the
///   options above. Requires `--taxon-ids`. When `--include-descendants` is set, each
///   excluded taxon's descendants are excluded as well. Excluded taxa that do not appear
///   in the report are ignored (with a warning), since a taxon with no reads in the
///   sample has nothing to exclude.
///
/// # Output
///
/// The number of `--output` (`-o`) paths sets the layout:
///
///   inputs   outputs   layout
///   1        1         single-end
///   2        2         R1 and R2 in separate files
///   2        1         interleaved (R1, R2, R1, R2, ...)
///
/// Each output's extension sets its compression: `.gz` and `.bgz` paths are written
/// bgzf-compressed, and any other path is written as uncompressed FASTQ. `-` writes
/// uncompressed FASTQ to stdout, and may be given for only one output.
///
/// `--threads` and `--compression-level` control bgzf compression and have no effect
/// when no output is bgzf-compressed.
///
/// # Examples
///
/// ```bash
/// # All reads classified as E. coli (taxon 562)
/// k2tools filter -r report.txt -k output.txt -i reads.fq.gz -o ecoli.fq.gz -t 562
///
/// # All Enterobacteriaceae (taxon 543), including every species and strain beneath it
/// k2tools filter -r report.txt -k output.txt -i reads.fq.gz -o entero.fq.gz -t 543 -d
///
/// # Unclassified reads from a paired-end run
/// k2tools filter -r report.txt -k output.txt \
///     -i r1.fq.gz r2.fq.gz -o unclass_r1.fq.gz unclass_r2.fq.gz -u
///
/// # E. coli read pairs, interleaved and uncompressed, streamed straight into an aligner
/// k2tools filter -r report.txt -k output.txt \
///     -i r1.fq.gz r2.fq.gz -o - -t 562 | bwa mem -p ref.fa -
///
/// # Human reads plus unclassified in a single pass
/// k2tools filter -r report.txt -k output.txt \
///     -i reads.fq.gz -o host_and_unclass.fq.gz -t 9606 -d -u
///
/// # All Felidae (taxon 9681) except the Panthera clade (taxon 9688)
/// k2tools filter -r report.txt -k output.txt \
///     -i reads.fq.gz -o cats.fq.gz -t 9681 -e 9688 -d
///
/// # Classified reads outside the human clade (taxon 9606). Reads assigned to an
/// # ancestor of human (e.g. Primates or Mammalia) are retained, so this is not a
/// # substitute for host depletion.
/// k2tools filter -r report.txt -k output.txt \
///     -i reads.fq.gz -o non_human.fq.gz -t 1 -e 9606 -d
/// ```
#[derive(clap::Args)]
#[command(verbatim_doc_comment)]
pub struct Filter {
    /// Path to the kraken2 report file.
    #[arg(short = 'r', long)]
    kraken_report: PathBuf,

    /// Path to the kraken2 per-read classification output.
    #[arg(short = 'k', long)]
    kraken_output: PathBuf,

    /// Input FASTQ file(s). One for single-end, two for paired-end.
    /// Supports gzip/bgzf compressed inputs.
    #[arg(short, long, num_args = 1..=2, required = true)]
    input: Vec<PathBuf>,

    /// Output FASTQ file(s): one per input, or one for paired-end input to interleave
    /// R1 and R2. `.gz`/`.bgz` paths are bgzf-compressed, other paths are uncompressed,
    /// and `-` writes uncompressed FASTQ to stdout.
    #[arg(short, long, num_args = 1..=2, required = true)]
    output: Vec<PathBuf>,

    /// Taxon ID(s) to extract reads for. At least one taxon ID or
    /// --include-unclassified must be specified.
    #[arg(short, long, num_args = 1..)]
    taxon_ids: Vec<u64>,

    /// Taxon ID(s) to exclude from the selected set. Requires --taxon-ids.
    /// Descendants are also excluded when --include-descendants is set.
    #[arg(short = 'e', long, num_args = 1..)]
    exclude_taxon_ids: Vec<u64>,

    /// Include reads assigned to any descendant of the specified taxa.
    #[arg(short = 'd', long, default_value_t = false)]
    include_descendants: bool,

    /// Include unclassified reads (taxon ID 0) in the output.
    #[arg(short = 'u', long, default_value_t = false)]
    include_unclassified: bool,

    /// Number of threads for bgzf compression. Unused if no output is bgzf.
    #[arg(long, default_value_t = 4)]
    threads: usize,

    /// Bgzf compression level (0-9). Unused if no output is bgzf.
    #[arg(long, default_value_t = 5)]
    compression_level: u8,
}

impl Command for Filter {
    fn execute(&self) -> Result<()> {
        self.validate_args()?;

        let report = KrakenReport::from_path(&self.kraken_report)?;
        if report.is_empty() {
            return self.handle_empty_inputs();
        }

        let (taxon_set, expected) = build_taxon_set_and_expected_count(
            &report,
            &self.taxon_ids,
            &self.exclude_taxon_ids,
            self.include_descendants,
            self.include_unclassified,
        )?;
        log::info!(
            "Filtering for {} taxa; expecting approximately {} reads",
            format_count(taxon_set.len() as u64),
            format_count(expected),
        );

        let (total, kept) = match self.run_filter_pipeline(&taxon_set) {
            Ok(counts) => counts,
            // A downstream reader closing stdout early (e.g. `| head`) is a normal way to
            // stop, not an input error
            Err(e) if self.output.iter().any(|p| is_stdout(p)) && is_broken_pipe(&e) => {
                log::info!("Stdout was closed by the downstream reader; stopping early.");
                return Ok(());
            }
            Err(e) => {
                self.print_error_banner(&e);
                return Err(e);
            }
        };

        #[allow(clippy::cast_precision_loss)]
        let pct = if total > 0 { kept as f64 / total as f64 * 100.0 } else { 0.0 };
        log::info!(
            "Kept {} / {} reads ({pct:.2}%), expected {}.",
            format_count(kept),
            format_count(total),
            format_count(expected),
        );

        Ok(())
    }
}

impl Filter {
    /// Handles the case where kraken2 was run on empty FASTQ files, producing an
    /// empty report and no kraken output file. Verifies that all FASTQ inputs are
    /// truly empty, then writes valid empty outputs.
    fn handle_empty_inputs(&self) -> Result<()> {
        let io = Io::new(u32::from(self.compression_level), IO_BUFFER_SIZE);
        for path in &self.input {
            let reader = io
                .new_reader(path)
                .with_context(|| format!("failed to open FASTQ: {}", path.display()))?;
            let mut fq = FastqReader::new(reader);
            if fq.next().is_some() {
                anyhow::bail!(
                    "kraken2 report is empty but FASTQ input {} contains records; \
                     inputs are inconsistent",
                    path.display()
                );
            }
        }

        let (pool, writers) = self.build_writers()?;
        close_writers(pool, writers)?;

        log::info!("Report is empty; all inputs are empty. Wrote empty outputs.");
        Ok(())
    }

    /// Validates command-line arguments beyond what clap enforces.
    fn validate_args(&self) -> Result<()> {
        let (num_inputs, num_outputs) = (self.input.len(), self.output.len());
        anyhow::ensure!(
            num_outputs == num_inputs || (num_inputs == 2 && num_outputs == 1),
            "got {num_inputs} input(s) and {num_outputs} output(s); give one output per input, \
             or a single output to interleave paired-end reads"
        );
        anyhow::ensure!(
            self.output.iter().filter(|p| is_stdout(p)).count() <= 1,
            "only one output may be written to stdout ('{STDOUT_PATH}')"
        );
        anyhow::ensure!(self.threads >= 1, "threads must be at least 1");
        anyhow::ensure!(self.compression_level <= 9, "compression level must be 0-9");
        anyhow::ensure!(
            !self.taxon_ids.is_empty() || self.include_unclassified,
            "at least one --taxon-ids value or --include-unclassified must be specified"
        );
        anyhow::ensure!(
            self.exclude_taxon_ids.is_empty() || !self.taxon_ids.is_empty(),
            "--exclude-taxon-ids requires --taxon-ids"
        );
        Ok(())
    }

    /// Opens all inputs, creates writers, runs the main filter loop, and closes
    /// everything down. Returns (total_reads, kept_reads).
    fn run_filter_pipeline(&self, taxon_set: &HashSet<u64>) -> Result<(u64, u64)> {
        let io = Io::new(u32::from(self.compression_level), IO_BUFFER_SIZE);
        let kraken_reader = io.new_reader(&self.kraken_output).with_context(|| {
            format!("failed to open kraken output: {}", self.kraken_output.display())
        })?;
        let mut kraken_iter = KrakenOutputReader::new(kraken_reader)
            .read_ahead(READ_AHEAD_CHUNK_SIZE, READ_AHEAD_NUM_CHUNKS);

        let is_paired = self.input.len() == 2;
        if is_paired && self.output.len() == 1 {
            log::info!("Writing paired-end reads interleaved to a single output");
        }
        let mut fq_iter1 = FastqReader::new(
            io.new_reader(&self.input[0])
                .with_context(|| format!("failed to open FASTQ: {}", self.input[0].display()))?,
        )
        .into_records()
        .read_ahead(READ_AHEAD_CHUNK_SIZE, READ_AHEAD_NUM_CHUNKS);

        let mut fq_iter2 = if is_paired {
            Some(
                FastqReader::new(io.new_reader(&self.input[1]).with_context(|| {
                    format!("failed to open FASTQ: {}", self.input[1].display())
                })?)
                .into_records()
                .read_ahead(READ_AHEAD_CHUNK_SIZE, READ_AHEAD_NUM_CHUNKS),
            )
        } else {
            None
        };

        let (pool, mut writers) = self.build_writers()?;
        let mut progress = ProgressLogger::new("k2tools::filter", "reads", 5_000_000);

        // Run the filter and verification, capturing any error so we can
        // shut down the pool cleanly before propagating it (avoids panics
        // in PooledWriter::drop when writers outlive the pool).
        let result = filter_reads(
            &mut kraken_iter,
            &mut fq_iter1,
            fq_iter2.as_mut(),
            taxon_set,
            &mut writers,
            &mut progress,
        )
        .and_then(|(total, kept)| {
            verify_fastq_exhausted(&mut fq_iter1, fq_iter2.as_mut(), total)?;
            Ok((total, kept))
        });

        progress.finish();

        // Close outputs even if filtering failed, but report the filtering error first
        let close_result = close_writers(pool, writers);
        let counts = result?;
        close_result?;
        Ok(counts)
    }

    /// Prints a prominent banner to stderr describing `error` and warning that the
    /// outputs may be incomplete.
    fn print_error_banner(&self, error: &anyhow::Error) {
        let banner = "#".repeat(72);
        let output_paths: Vec<_> =
            self.output
                .iter()
                .map(|p| {
                    if is_stdout(p) { "  stdout".to_string() } else { format!("  {}", p.display()) }
                })
                .collect();
        eprintln!(
            "\n{banner}\n\
             # ERROR: invalid inputs detected\n\
             #\n\
             # {error}\n\
             #\n\
             # WARNING: partial/invalid output files may have been written to:\n\
             # {}\n\
             {banner}\n",
            output_paths.join("\n"),
        );
    }

    /// Opens a writer for each output path. `-` writes uncompressed to stdout, `.gz` and
    /// `.bgz` paths are bgzf-compressed through a writer pool, and any other path is
    /// written uncompressed. The pool is only created if at least one output is bgzf.
    ///
    /// Returns (pool, writers) so that destructuring as `let (pool, writers) = ...`
    /// ensures writers are dropped before the pool (reverse declaration order).
    fn build_writers(&self) -> Result<(Option<Pool>, Vec<FastqWriter>)> {
        let mut pool_builder = if self.output.iter().any(Io::is_gzip_path) {
            Some(
                PoolBuilder::<BufWriter<File>, BgzfCompressor>::new()
                    .threads(self.threads)
                    .queue_size(self.threads * 50)
                    .compression_level(self.compression_level)?,
            )
        } else {
            None
        };

        let mut writers = Vec::with_capacity(self.output.len());
        for path in &self.output {
            if is_stdout(path) {
                writers.push(FastqWriter::plain(Box::new(std::io::stdout().lock())));
                continue;
            }
            let file = File::create(path)
                .with_context(|| format!("failed to create output: {}", path.display()))?;
            let writer = match pool_builder.as_mut() {
                Some(builder) if Io::is_gzip_path(path) => {
                    FastqWriter::Bgzf(builder.exchange(BufWriter::new(file)))
                }
                _ => FastqWriter::plain(Box::new(file)),
            };
            writers.push(writer);
        }

        let pool = pool_builder.map(PoolBuilder::build).transpose()?;
        Ok((pool, writers))
    }
}

/// A FASTQ output: bgzf-compressed through the shared writer pool, or uncompressed and
/// written directly on the calling thread.
enum FastqWriter {
    Bgzf(PooledWriter),
    Plain(BufWriter<Box<dyn Write>>),
}

impl FastqWriter {
    /// Wraps `inner` in a buffered, uncompressed writer.
    fn plain(inner: Box<dyn Write>) -> Self {
        FastqWriter::Plain(BufWriter::with_capacity(IO_BUFFER_SIZE, inner))
    }

    /// Flushes any buffered data and closes the writer.
    fn close(self) -> Result<()> {
        match self {
            FastqWriter::Bgzf(writer) => writer.close()?,
            FastqWriter::Plain(mut writer) => writer.flush()?,
        }
        Ok(())
    }
}

impl Write for FastqWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            FastqWriter::Bgzf(writer) => writer.write(buf),
            FastqWriter::Plain(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            FastqWriter::Bgzf(writer) => writer.flush(),
            FastqWriter::Plain(writer) => writer.flush(),
        }
    }
}

/// Closes all writers and then stops the pool, if there is one. Pooled writers must be
/// closed before their pool is stopped. Every writer is closed and the pool stopped even
/// if an earlier step fails, and the first error is returned.
fn close_writers(pool: Option<Pool>, writers: Vec<FastqWriter>) -> Result<()> {
    let mut first_error = None;
    for writer in writers {
        if let Err(e) = writer.close() {
            first_error.get_or_insert(e);
        }
    }
    if let Some(mut pool) = pool {
        if let Err(e) = pool.stop_pool() {
            first_error.get_or_insert(e.into());
        }
    }
    first_error.map_or(Ok(()), Err)
}

/// Returns true if `error` was caused by writing to a pipe whose reader has closed.
fn is_broken_pipe(error: &anyhow::Error) -> bool {
    error.chain().any(|cause| {
        cause
            .downcast_ref::<std::io::Error>()
            .is_some_and(|e| e.kind() == std::io::ErrorKind::BrokenPipe)
    })
}

/// Returns true if `path` is the stdout marker `-`.
fn is_stdout(path: &Path) -> bool {
    path == Path::new(STDOUT_PATH)
}

/// Runs the main filter loop: co-iterates kraken output and FASTQ iterator(s) in
/// lockstep, writing matching records to the output writers. R2 records go to the
/// second writer if there are two, or are interleaved after R1 if there is only one.
///
/// Returns (total_reads_processed, reads_kept).
fn filter_reads(
    kraken_iter: &mut impl Iterator<Item = Result<KrakenRecord>>,
    fq_iter1: &mut impl Iterator<Item = Result<OwnedRecord, FastqError>>,
    mut fq_iter2: Option<&mut impl Iterator<Item = Result<OwnedRecord, FastqError>>>,
    taxon_set: &HashSet<u64>,
    writers: &mut [FastqWriter],
    progress: &mut ProgressLogger,
) -> Result<(u64, u64)> {
    let mut total: u64 = 0;
    let mut kept: u64 = 0;
    let r2_writer_index = writers.len() - 1;

    for kraken_result in kraken_iter {
        let kraken_rec = kraken_result?;
        total += 1;
        progress.record();

        let fq_rec1 = fq_iter1
            .next()
            .context("FASTQ input ended before kraken output")?
            .with_context(|| format!("failed to read FASTQ record at kraken line {total}"))?;

        let fq_rec2: Option<OwnedRecord> = if let Some(ref mut iter2) = fq_iter2 {
            Some(
                iter2
                    .next()
                    .context("second FASTQ input ended before kraken output")?
                    .with_context(|| {
                        format!("failed to read FASTQ R2 record at kraken line {total}")
                    })?,
            )
        } else {
            None
        };

        if taxon_set.contains(&kraken_rec.taxon_id()) {
            // Validate read names only for matching reads to avoid overhead
            validate_read_name(kraken_rec.read_name(), fq_rec1.head(), total)?;
            if let Some(ref rec2) = fq_rec2 {
                validate_read_name(kraken_rec.read_name(), rec2.head(), total)?;
            }

            write_fastq_record(&mut writers[0], &fq_rec1)?;
            if let Some(ref rec2) = fq_rec2 {
                write_fastq_record(&mut writers[r2_writer_index], rec2)?;
            }
            kept += 1;
        }
    }

    Ok((total, kept))
}

/// Verifies that the FASTQ streams are exhausted after the kraken output ends.
fn verify_fastq_exhausted(
    fq_iter1: &mut impl Iterator<Item = Result<OwnedRecord, FastqError>>,
    fq_iter2: Option<&mut impl Iterator<Item = Result<OwnedRecord, FastqError>>>,
    total: u64,
) -> Result<()> {
    if fq_iter1.next().is_some() {
        anyhow::bail!("FASTQ input has more records than kraken output ({total} kraken records)");
    }
    if let Some(iter2) = fq_iter2 {
        if iter2.next().is_some() {
            anyhow::bail!(
                "second FASTQ input has more records than kraken output ({total} kraken records)"
            );
        }
    }
    Ok(())
}

/// Builds the set of taxon IDs to filter for and computes the expected number of
/// matching reads from the report's count fields.
///
/// If `include_descendants` is true, expands each taxon ID (included and excluded)
/// to cover all its descendants in the report taxonomy tree. If `include_unclassified`
/// is true, adds taxon ID 0. Taxa in `exclude_taxon_ids` are then removed from the
/// set; excluded taxa not present in the report are ignored with a warning, since a
/// taxon with no reads in the sample has nothing to exclude, and an exclusion that
/// removes nothing from the selection warns (usually a forgotten `-d`). The expected
/// count is the sum of `direct_count` over the final set.
///
/// Returns `(taxon_id_set, expected_read_count)`. Errors if exclusion removes every
/// selected taxon, since that would silently produce empty outputs.
fn build_taxon_set_and_expected_count(
    report: &KrakenReport,
    taxon_ids: &[u64],
    exclude_taxon_ids: &[u64],
    include_descendants: bool,
    include_unclassified: bool,
) -> Result<(HashSet<u64>, u64)> {
    let mut set = HashSet::new();

    for &tid in taxon_ids {
        let idx = report
            .index_of_taxon_id(tid)
            .with_context(|| format!("taxon ID {tid} not found in report"))?;
        set.insert(tid);
        if include_descendants {
            for desc_idx in report.descendants(idx) {
                set.insert(report.row(desc_idx).taxon_id());
            }
        }
    }

    if include_unclassified {
        set.insert(0);
    }

    for &tid in exclude_taxon_ids {
        let Some(idx) = report.index_of_taxon_id(tid) else {
            log::warn!("Excluded taxon ID {tid} not found in report; ignoring");
            continue;
        };

        let mut indices = vec![idx];
        if include_descendants {
            indices.extend(report.descendants(idx));
        }
        let mut removed_any = false;
        for i in indices {
            removed_any |= set.remove(&report.row(i).taxon_id());
        }
        if !removed_any {
            log::warn!(
                "Excluded taxon ID {tid} removed nothing from the selection; it was not \
                 selected by --taxon-ids (missing --include-descendants?) or already excluded"
            );
        }
    }

    anyhow::ensure!(
        !set.is_empty(),
        "--exclude-taxon-ids removed every selected taxon; nothing would be extracted"
    );

    // Each read is counted in the direct count of exactly one taxon, so this is exact
    // even when selected clades overlap (e.g. `-t 2 543 -d`)
    let expected = set
        .iter()
        .filter_map(|&tid| report.get_by_taxon_id(tid))
        .map(ReportRow::direct_count)
        .sum();

    Ok((set, expected))
}

/// Validates that a kraken read name matches a FASTQ record header.
///
/// Expects the FASTQ header to start with the kraken read name (byte-for-byte),
/// optionally followed by `/1` or `/2` (paired-end suffix) and/or whitespace
/// plus a comment. Avoids scanning the full header — only checks the prefix
/// at the kraken name length boundary.
fn validate_read_name(kraken_name: &str, fastq_head: &[u8], line_number: u64) -> Result<()> {
    let k = kraken_name.as_bytes();
    let f = fastq_head;

    if f.len() >= k.len() && f[..k.len()] == *k {
        let rest = &f[k.len()..];
        if rest.is_empty()
            || rest[0] == b' '
            || rest[0] == b'\t'
            || (rest.len() >= 2
                && rest[0] == b'/'
                && (rest[1] == b'1' || rest[1] == b'2')
                && (rest.len() == 2 || rest[2] == b' ' || rest[2] == b'\t'))
        {
            return Ok(());
        }
    }

    // Build a readable FASTQ name for the error message only on failure
    let name_end = f.iter().position(|&b| b == b' ' || b == b'\t').unwrap_or(f.len());
    anyhow::bail!(
        "read name mismatch at kraken line {line_number}: \
         kraken={kraken_name:?}, FASTQ={:?}",
        String::from_utf8_lossy(&f[..name_end])
    );
}

/// Writes a single FASTQ record to a writer.
fn write_fastq_record<W: Write>(writer: &mut W, rec: &impl Record) -> Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(rec.head())?;
    writer.write_all(b"\n")?;
    writer.write_all(rec.seq())?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(rec.qual())?;
    writer.write_all(b"\n")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_report() -> KrakenReport {
        // unclassified(0), root(1), Bacteria(2), E.coli(3), Eukaryota(4), Human(5).
        // Counts are internally consistent: each clade count equals the taxon's
        // direct count plus its children's clade counts, as kraken2 guarantees.
        let lines = [
            " 10.00\t100\t100\tU\t0\tunclassified",
            " 90.00\t900\t0\tR\t1\troot",
            " 60.00\t600\t100\tD\t2\t  Bacteria",
            " 50.00\t500\t500\tS\t3\t    Escherichia coli",
            " 30.00\t300\t100\tD\t4\t  Eukaryota",
            " 20.00\t200\t200\tS\t5\t    Homo sapiens",
        ]
        .join("\n");
        KrakenReport::from_reader(lines.as_bytes()).unwrap()
    }

    #[test]
    fn test_build_taxon_set_exact() {
        let report = make_report();
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[3], &[], false, false).unwrap();
        assert_eq!(set, HashSet::from([3]));
        assert_eq!(expected, 500);
    }

    #[test]
    fn test_build_taxon_set_with_descendants() {
        let report = make_report();
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[2], &[], true, false).unwrap();
        assert_eq!(set, HashSet::from([2, 3]));
        assert_eq!(expected, 600);
    }

    #[test]
    fn test_build_taxon_set_with_descendants_root() {
        let report = make_report();
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[1], &[], true, false).unwrap();
        assert_eq!(set, HashSet::from([1, 2, 3, 4, 5]));
        assert_eq!(expected, 900);
    }

    #[test]
    fn test_build_taxon_set_unknown_taxon() {
        let report = make_report();
        let result = build_taxon_set_and_expected_count(&report, &[99999], &[], false, false);
        assert!(result.is_err());
    }

    #[test]
    fn test_build_taxon_set_include_unclassified() {
        let report = make_report();
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[3], &[], false, true).unwrap();
        assert_eq!(set, HashSet::from([0, 3]));
        assert_eq!(expected, 600);
    }

    #[test]
    fn test_build_taxon_set_only_unclassified() {
        let report = make_report();
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[], &[], false, true).unwrap();
        assert_eq!(set, HashSet::from([0]));
        assert_eq!(expected, 100);
    }

    #[test]
    fn test_expected_count_with_descendants() {
        let report = make_report();
        let (_, expected) =
            build_taxon_set_and_expected_count(&report, &[2], &[], true, false).unwrap();
        assert_eq!(expected, 600);
    }

    #[test]
    fn test_expected_count_without_descendants() {
        let report = make_report();
        let (_, expected) =
            build_taxon_set_and_expected_count(&report, &[2], &[], false, false).unwrap();
        assert_eq!(expected, 100);
    }

    #[test]
    fn test_expected_count_with_unclassified() {
        let report = make_report();
        let (_, expected) =
            build_taxon_set_and_expected_count(&report, &[3], &[], false, true).unwrap();
        assert_eq!(expected, 600);
    }

    #[test]
    fn test_expected_count_does_not_double_count_overlapping_clades() {
        let report = make_report();
        // E.coli (3) is inside the Bacteria (2) clade, so its reads must count once
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[2, 3], &[], true, false).unwrap();
        assert_eq!(set, HashSet::from([2, 3]));
        assert_eq!(expected, 600);
    }

    #[test]
    fn test_exclude_removes_taxon_and_its_expected_count() {
        let report = make_report();
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[3, 5], &[5], false, false).unwrap();
        assert_eq!(set, HashSet::from([3]));
        assert_eq!(expected, 500);
    }

    #[test]
    fn test_exclude_with_descendants_removes_whole_clade() {
        let report = make_report();
        // Everything under root except the Eukaryota clade (4 and Human 5); the
        // expected count is exactly clade(root) - clade(Eukaryota)
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[1], &[4], true, false).unwrap();
        assert_eq!(set, HashSet::from([1, 2, 3]));
        assert_eq!(expected, 900 - 300);
    }

    #[test]
    fn test_exclude_of_taxon_not_in_set_is_noop() {
        let report = make_report();
        // E.coli (3) is a descendant of Bacteria (2) but descendants were not included
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[2], &[3], false, false).unwrap();
        assert_eq!(set, HashSet::from([2]));
        assert_eq!(expected, 100);
    }

    #[test]
    fn test_exclude_taxon_missing_from_report_is_ignored() {
        let report = make_report();
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[3], &[99999], false, false).unwrap();
        assert_eq!(set, HashSet::from([3]));
        assert_eq!(expected, 500);
    }

    #[test]
    fn test_exclude_of_other_taxon_leaves_unclassified_in_set() {
        let report = make_report();
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[1], &[5], true, true).unwrap();
        assert_eq!(set, HashSet::from([0, 1, 2, 3, 4]));
        assert_eq!(expected, 900 + 100 - 200);
    }

    #[test]
    fn test_exclude_of_taxon_zero_removes_unclassified() {
        let report = make_report();
        let (set, expected) =
            build_taxon_set_and_expected_count(&report, &[1], &[0], true, true).unwrap();
        assert_eq!(set, HashSet::from([1, 2, 3, 4, 5]));
        assert_eq!(expected, 900);
    }

    #[test]
    fn test_exclusion_of_all_selected_taxa_errors() {
        let report = make_report();
        let result = build_taxon_set_and_expected_count(&report, &[3], &[3], false, false);
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_read_name_match() {
        assert!(validate_read_name("read1", b"read1", 1).is_ok());
    }

    #[test]
    fn test_validate_read_name_mismatch() {
        assert!(validate_read_name("read1", b"read2", 1).is_err());
    }

    #[test]
    fn test_validate_read_name_strip_suffix_1() {
        assert!(validate_read_name("read1", b"read1/1", 1).is_ok());
    }

    #[test]
    fn test_validate_read_name_strip_suffix_2() {
        assert!(validate_read_name("read1", b"read1/2", 1).is_ok());
    }

    #[test]
    fn test_validate_read_name_with_comment() {
        assert!(validate_read_name("read1", b"read1 length=150", 1).is_ok());
    }

    #[test]
    fn test_validate_read_name_suffix_and_comment() {
        assert!(validate_read_name("read1", b"read1/1 length=150", 1).is_ok());
    }

    /// Builds an otherwise-valid `Filter` with the given input and output paths.
    fn filter_with_io(inputs: &[&str], outputs: &[&str]) -> Filter {
        Filter {
            kraken_report: PathBuf::from("r.txt"),
            kraken_output: PathBuf::from("k.txt"),
            input: inputs.iter().map(PathBuf::from).collect(),
            output: outputs.iter().map(PathBuf::from).collect(),
            taxon_ids: vec![1],
            exclude_taxon_ids: vec![],
            include_descendants: false,
            include_unclassified: false,
            threads: 4,
            compression_level: 6,
        }
    }

    #[test]
    fn test_validate_args_rejects_two_outputs_for_single_input() {
        let filter = filter_with_io(&["a.fq"], &["b.fq", "c.fq"]);
        assert!(filter.validate_args().is_err());
    }

    #[test]
    fn test_validate_args_allows_single_output_for_paired_input() {
        let filter = filter_with_io(&["a.fq", "b.fq"], &["c.fq"]);
        assert!(filter.validate_args().is_ok());
    }

    #[test]
    fn test_is_broken_pipe_detects_broken_pipe_io_error() {
        let error = anyhow::Error::from(std::io::Error::from(std::io::ErrorKind::BrokenPipe));
        assert!(is_broken_pipe(&error));
    }

    #[test]
    fn test_is_broken_pipe_detects_broken_pipe_under_context() {
        let error = anyhow::Error::from(std::io::Error::from(std::io::ErrorKind::BrokenPipe))
            .context("failed to write output");
        assert!(is_broken_pipe(&error));
    }

    #[test]
    fn test_is_broken_pipe_ignores_other_io_errors() {
        let error = anyhow::Error::from(std::io::Error::from(std::io::ErrorKind::NotFound));
        assert!(!is_broken_pipe(&error));
    }

    #[test]
    fn test_is_broken_pipe_ignores_non_io_errors() {
        assert!(!is_broken_pipe(&anyhow::anyhow!("read name mismatch")));
    }

    #[test]
    fn test_validate_args_rejects_stdout_for_both_outputs() {
        let filter = filter_with_io(&["a.fq", "b.fq"], &["-", "-"]);
        assert!(filter.validate_args().is_err());
    }

    #[test]
    fn test_validate_args_no_taxa_or_unclassified() {
        let filter = Filter {
            kraken_report: PathBuf::from("r.txt"),
            kraken_output: PathBuf::from("k.txt"),
            input: vec![PathBuf::from("a.fq")],
            output: vec![PathBuf::from("b.fq")],
            taxon_ids: vec![],
            exclude_taxon_ids: vec![],
            include_descendants: false,
            include_unclassified: false,
            threads: 4,
            compression_level: 6,
        };
        assert!(filter.validate_args().is_err());
    }

    #[test]
    fn test_validate_args_exclude_requires_taxon_ids() {
        let filter = Filter {
            kraken_report: PathBuf::from("r.txt"),
            kraken_output: PathBuf::from("k.txt"),
            input: vec![PathBuf::from("a.fq")],
            output: vec![PathBuf::from("b.fq")],
            taxon_ids: vec![],
            exclude_taxon_ids: vec![5],
            include_descendants: false,
            include_unclassified: true,
            threads: 4,
            compression_level: 6,
        };
        assert!(filter.validate_args().is_err());
    }
}
