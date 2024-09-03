#[cfg(test)]
mod test;

mod ode1be;
mod ode23;
mod ode45;

pub use ode1be::ode1be;
pub use ode23::ode23;
pub use ode45::ode45;

use crate::get_defeat_message;
use std::fmt;

/// Possible errors encountered when integrating.
pub enum IntegrationError {
    GeneralError,
}

/// Debug implementation for solve errors.
impl fmt::Debug for IntegrationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let error = match self {
            Self::GeneralError => "\x1b[1;91m???.\x1b[0;91m".to_string(),
        };
        write!(
            f,
            "\n{}\n\x1b[0;2;31m{}\x1b[0m\n",
            error,
            get_defeat_message()
        )
    }
}

/// Display implementation for solve errors.
impl fmt::Display for IntegrationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let error = match self {
            Self::GeneralError => "\x1b[1;91m???.\x1b[0;91m".to_string(),
        };
        write!(
            f,
            "\n{}\n\x1b[0;2;31m{}\x1b[0m\n",
            error,
            get_defeat_message()
        )
    }
}

/// Implementation of solve errors from ok_or(&str).
impl From<&str> for IntegrationError {
    fn from(string: &str) -> Self {
        todo!("{}", string)
    }
}
