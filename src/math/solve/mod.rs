#[cfg(test)]
mod test;

mod ode23;

use crate::get_defeat_message;
use std::fmt;

/// ???
pub enum SolveError {
    GeneralError,
}

/// ???
impl fmt::Debug for SolveError {
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

/// ???
impl fmt::Display for SolveError {
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

/// ???
impl From<&str> for SolveError {
    fn from(string: &str) -> Self {
        todo!("{}", string)
    }
}
