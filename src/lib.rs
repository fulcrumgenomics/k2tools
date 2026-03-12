#![deny(unsafe_code)]

#[cfg(any(target_pointer_width = "16", target_pointer_width = "32"))]
compile_error!("k2tools requires a 64-bit or wider platform");

pub mod commands;
pub mod kraken_output;
pub mod progress;
pub mod report;
pub mod version;
