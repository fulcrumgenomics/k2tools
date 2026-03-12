# Contributing to k2tools

k2tools is a Rust CLI toolkit for working with
[kraken2](https://github.com/DerrickWood/kraken2) outputs. It provides
commands for filtering, summarizing, and manipulating kraken2 classification
results and report files.

## Getting Started

**Prerequisites:** Rust stable (see `Cargo.toml` for minimum version) and
[cargo-nextest](https://nexte.st/).

```bash
cargo build              # debug build
cargo build --release    # optimized build
```

### Verification Checklist

Run all three before submitting changes:

```bash
cargo ci-fmt             # check formatting
cargo ci-lint            # clippy, pedantic, warnings-as-errors
cargo ci-test            # all tests via nextest
```

The `ci-*` aliases are defined in `.cargo/config.toml`.

## Architecture

### Commands

Each subcommand is a struct that implements the `Command` trait
(`src/commands/command.rs`) and is dispatched via the `Subcommand` enum in
`src/main.rs`.

## Adding a New Command

1. Create `src/commands/<name>.rs` containing a command struct (derives
   `clap::Args`, implements `Command`)
2. Add `pub mod <name>;` to `src/commands/mod.rs`
3. Add a variant to `Subcommand` in `src/main.rs` and a match arm in `execute()`
4. Add tests
5. Run the full verification checklist

## Code Style

**Priorities:** correctness > readability > performance.

- Write idiomatic Rust — don't fight the language
- Prefer meaningful names even if longer (`taxon_id_set` over `set`)
- Extract small/medium functions with clear inputs and outputs
- Keep related code together
- Use `pub(crate)` for shared internal utilities
- Doc comments on all public items; code comments explain *why*, not *what*
- Don't over-generalize — solve the problem at hand

## Testing

**No checked-in test data.** All test data is built programmatically.

**Guidelines:**

- Integration tests go in `tests/test_<command>.rs`; unit tests in inline
  `#[cfg(test)]` modules
- Prefer many small, focused tests over parameterized/table-driven tests
- Test function, not implementation — tests should survive a significant refactor
- Cover expected results, error conditions, and boundary cases

## Output Format Conventions

- **TSV output**: lowercase `snake_case` headers, tab-separated, no metadata or
  comment lines
- **Fractions**: use `frac_` prefix (not `pct_`)

## Reporting Bugs

Please include:

1. **Version**: output of `k2tools --version`
2. **Command**: the exact command line that failed
3. **Error output**: full stderr/stdout
4. **OS**: operating system and version
5. **Expected vs. actual behavior**: what should happen vs. what did happen
