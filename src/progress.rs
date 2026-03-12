use std::time::Instant;

use log::Level;

/// Width of the formatted count field, enough for `"1,000,000,000"` (13 chars).
const COUNT_WIDTH: usize = 13;

/// Format a count with commas as thousands separators (e.g. `1000000` → `"1,000,000"`).
pub(crate) fn format_count(n: u64) -> String {
    let s = n.to_string();
    s.as_bytes()
        .rchunks(3)
        .rev()
        .map(|c| std::str::from_utf8(c).unwrap())
        .collect::<Vec<_>>()
        .join(",")
}

/// Format a duration as `MMm SSs` with zero-padded two-digit minutes and seconds.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "elapsed seconds are non-negative and fit in u64 for any practical duration"
)]
fn format_elapsed(secs: f64) -> String {
    let total_secs = secs as u64;
    let m = total_secs / 60;
    let s = total_secs % 60;
    format!("{m:02}m {s:02}s")
}

/// Logs progress every N records.
///
/// Designed for high-throughput loops where millions of items are processed.
/// The counter is checked against a pre-computed milestone to avoid a division
/// on every call; only when the milestone is reached does the logger perform
/// any formatting or I/O.
pub struct ProgressLogger {
    name: &'static str,
    unit: &'static str,
    every: u64,
    count: u64,
    /// Next count at which to emit a progress message. Avoids per-call division.
    next_milestone: u64,
    last_milestone: Instant,
    start: Instant,
}

impl ProgressLogger {
    /// Creates a new progress logger.
    ///
    /// * `name` — log target (appears in `[target]` in log output)
    /// * `unit` — label for items being counted (e.g. `"reads"`, `"records"`)
    /// * `every` — emit a progress line every N items
    #[must_use]
    pub fn new(name: &'static str, unit: &'static str, every: u64) -> Self {
        let now = Instant::now();
        Self { name, unit, every, count: 0, next_milestone: every, last_milestone: now, start: now }
    }

    /// Record one item and log if the interval has been reached.
    pub fn record(&mut self) {
        self.count += 1;
        if self.count >= self.next_milestone {
            self.next_milestone += self.every;
            self.emit();
        }
    }

    /// Advance the counter by `n` and emit a progress message if a milestone
    /// boundary is crossed. If `n` spans multiple milestones only one message
    /// is emitted.
    pub fn record_n(&mut self, n: u64) {
        if n == 0 {
            return;
        }
        self.count += n;
        if self.count >= self.next_milestone {
            while self.next_milestone <= self.count {
                self.next_milestone += self.every;
            }
            self.emit();
        }
    }

    /// Log final totals.
    pub fn finish(&self) {
        let total = format_elapsed(self.start.elapsed().as_secs_f64());
        log::log!(
            target: self.name, Level::Info,
            "Processed {:>COUNT_WIDTH$} {} total in {total}.",
            format_count(self.count), self.unit,
        );
    }

    /// Emit a progress line with timing information.
    fn emit(&mut self) {
        let milestone_secs = self.last_milestone.elapsed().as_secs_f64();
        let total_elapsed = format_elapsed(self.start.elapsed().as_secs_f64());
        let last_took = format!("last {} took {:.1}s", format_count(self.every), milestone_secs);

        log::log!(
            target: self.name, Level::Info,
            "Processed {:>COUNT_WIDTH$} {} - elapsed time {total_elapsed} - {last_took}.",
            format_count(self.count), self.unit,
        );
        self.last_milestone = Instant::now();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_count_small() {
        assert_eq!(format_count(0), "0");
        assert_eq!(format_count(1), "1");
        assert_eq!(format_count(999), "999");
    }

    #[test]
    fn test_format_count_thousands() {
        assert_eq!(format_count(1_000), "1,000");
        assert_eq!(format_count(10_000), "10,000");
        assert_eq!(format_count(100_000), "100,000");
    }

    #[test]
    fn test_format_count_millions() {
        assert_eq!(format_count(1_000_000), "1,000,000");
        assert_eq!(format_count(10_000_000), "10,000,000");
        assert_eq!(format_count(1_234_567_890), "1,234,567,890");
    }

    #[test]
    fn test_format_elapsed_seconds_only() {
        assert_eq!(format_elapsed(5.3), "00m 05s");
    }

    #[test]
    fn test_format_elapsed_minutes_and_seconds() {
        assert_eq!(format_elapsed(125.7), "02m 05s");
    }

    #[test]
    fn test_format_elapsed_exact_minute() {
        assert_eq!(format_elapsed(60.0), "01m 00s");
    }

    #[test]
    fn test_record_no_panic() {
        let mut pl = ProgressLogger::new("test", "reads", 10);
        for _ in 0..25 {
            pl.record();
        }
        pl.finish();
    }

    #[test]
    fn test_record_n_zero_is_noop() {
        let mut pl = ProgressLogger::new("test", "reads", 10);
        pl.record_n(0);
        pl.record_n(0);
        pl.record_n(0);
        pl.finish();
    }

    #[test]
    fn test_record_n_crosses_milestones() {
        let mut pl = ProgressLogger::new("test", "reads", 10);
        pl.record_n(35);
        pl.finish();
    }
}
