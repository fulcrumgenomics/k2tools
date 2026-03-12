use anyhow::Result;

/// Trait implemented by every k2tools subcommand.
pub trait Command {
    /// Execute the command.
    ///
    /// # Errors
    /// Returns an error if the command fails.
    fn execute(&self) -> Result<()>;
}
