#![deny(unsafe_code)]

use std::io::Write;
use std::time::Instant;

use anyhow::Result;
use chrono::Local;
use clap::Parser;
use clap::builder::styling::{AnsiColor, Effects, Styles};
use env_logger::Env;
use k2tools_lib::commands::command::Command;
use k2tools_lib::commands::filter::Filter;
use k2tools_lib::commands::report_to_tsv::ReportToTsv;

#[global_allocator]
static GLOBAL: rpmalloc::RpMalloc = rpmalloc::RpMalloc;

const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Cyan.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Cyan.on_default())
    .valid(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .invalid(AnsiColor::Yellow.on_default().effects(Effects::BOLD));

/// Toolkit for working with kraken2 outputs.
///
/// k2tools provides commands for filtering, summarizing, and
/// manipulating the output files produced by kraken2.
#[derive(Parser)]
#[command(
    name = "k2tools",
    version = env!("CARGO_PKG_VERSION"),
    long_about,
    styles = STYLES,
    after_long_help = "Run 'k2tools <command> --help' for detailed options on each command."
)]
struct Cli {
    /// Enable verbose logging.
    #[arg(short, long, global = true)]
    verbose: bool,

    #[command(subcommand)]
    command: Subcommand,
}

/// All k2tools subcommands.
#[derive(clap::Subcommand)]
enum Subcommand {
    Filter(Filter),
    ReportToTsv(ReportToTsv),
}

impl Command for Subcommand {
    fn execute(&self) -> Result<()> {
        match self {
            Subcommand::Filter(c) => c.execute(),
            Subcommand::ReportToTsv(c) => c.execute(),
        }
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let log_level = if cli.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(Env::default().default_filter_or(log_level))
        .format(|buf, record| {
            let timestamp = Local::now().format("%Y-%m-%d %H:%M:%S");
            writeln!(buf, "[{timestamp} {} {}] {}", record.level(), record.target(), record.args())
        })
        .init();

    let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
    log::info!("k2tools by Fulcrum Genomics | https://https://github.com/fulcrumgenomics/k2tools");
    log::info!("Executing: {cmdline}");

    let start = Instant::now();
    let result = cli.command.execute();
    let elapsed = start.elapsed();
    let minutes = elapsed.as_secs() / 60;
    let seconds = elapsed.as_secs() % 60;

    match &result {
        Ok(()) => log::info!("Successfully completed execution in {minutes}m:{seconds:02}s."),
        Err(_) => log::error!("Execution failed after {minutes}m:{seconds:02}s."),
    }

    result
}
