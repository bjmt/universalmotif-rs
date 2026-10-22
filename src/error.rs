use thiserror::Error;

/// Errors returned by this crate.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum Error {
    /// A motif failed validation; each entry describes one problem.
    #[error("invalid motif:\n{}", .0.join("\n"))]
    InvalidMotif(Vec<String>),

    /// A function argument was unusable.
    #[error("{0}")]
    Argument(String),

    /// The requested operation is not supported for this input.
    #[error("{0}")]
    Unsupported(String),
}

pub type Result<T> = std::result::Result<T, Error>;

pub(crate) fn arg_err<T>(msg: impl Into<String>) -> Result<T> {
    Err(Error::Argument(msg.into()))
}
