use super::super::TensorRank0;
use crate::{ABS_TOL, REL_TOL};

/// ???
#[derive(Debug)]
pub struct Newton {
    /// Absolute error tolerance.
    pub _abs_tol: TensorRank0,
    /// Relative error tolerance.
    pub _rel_tol: TensorRank0,
}

impl Default for Newton {
    fn default() -> Self {
        Self {
            _abs_tol: ABS_TOL,
            _rel_tol: REL_TOL,
        }
    }
}

// scipy has an rtol, how it it used?
// and it doesnt have an atol, but think we can use it still
