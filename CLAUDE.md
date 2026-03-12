# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

k2tools is a Rust CLI toolkit for working with [kraken2](https://github.com/DerrickWood/kraken2) outputs. It provides commands for filtering, summarizing, and manipulating the classification results and report files produced by kraken2.

## Commands

```bash
# Build
cargo build
cargo build --release

# Test (uses nextest)
cargo ci-test                          # all tests
cargo nextest run <test_name>          # single test by name

# Lint and format
cargo ci-lint                          # clippy with pedantic warnings as errors
cargo ci-fmt                           # check formatting (--check mode)
cargo fmt --package k2tools            # apply formatting
```

The `ci-*` aliases are defined in `.cargo/config.toml`. CI runs all three.

## Architecture

### Crate Structure

- `k2tools` — binary + library (`k2tools_lib`)

The library is named `k2tools_lib` (see `[lib]` in `Cargo.toml`). External code and tests reference it as `k2tools_lib::...`.

### Command Pattern

Each subcommand is a struct implementing the `Command` trait (`src/commands/command.rs`). The `Subcommand` enum in `main.rs` dispatches to them via a `match`.

To add a new command:
1. Create `src/commands/<name>.rs` with the command struct
2. Add `pub mod <name>;` to `src/commands/mod.rs`
3. Add a variant to `Subcommand` in `src/main.rs` and a match arm in `Command::execute()`
4. Add tests

## Testing Conventions

- **No checked-in test data** — all test data is built programmatically
- Integration tests live in `tests/test_<command>.rs`; unit tests are inline `#[cfg(test)]` modules
- Prefer many small, focused tests over parameterized/table-driven tests
- Test function, not implementation — tests should survive a significant refactor

## Output Format Conventions

- TSV files: lowercase snake_case headers, tab-separated, no metadata or comment lines
- Fractions use `frac_` prefix (not `pct_`)

## Code Style

**Priorities:** correctness > readability > performance.

- Write idiomatic Rust
- Prefer meaningful names even if longer
- Extract small/medium functions with clear inputs and outputs
- Keep related code together
- Use `pub(crate)` for shared internal utilities
- Doc comments on all public items; code comments explain *why*, not *what*
- Don't over-generalize — solve the problem at hand
