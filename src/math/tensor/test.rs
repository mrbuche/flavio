use super::Tensor;
use crate::{get_defeat_message, ABS_TOL, REL_TOL};
use std::{cmp::PartialEq, fmt};

pub fn assert_eq<'a, T: PartialEq + Tensor>(
    value_1: &'a T,
    value_2: &'a T,
) -> Result<(), TestError> {
    if value_1 == value_2 {
        Ok(())
    } else {
        Err(TestError {
            message: format!(
            "\n\x1b[1;91mAssertion `left == right` failed.\n\x1b[0;91m  left: {}\n right: {}\x1b[0",
            value_1, value_2
        ),
        })
    }
}

pub fn assert_eq_within_tols<'a, T: Tensor>(
    value_1: &'a T,
    value_2: &'a T,
) -> Result<(), TestError> {
    //
    // Rank-4 panicking because norm_squared() called by norm() panics.
    //
    // I feel like this is more fair to compare entry-by-entry.
    //
    // So maybe you can require something like iter_flat() for Tensor,
    // and then use it here to do entry-by-entry rank-independent checks.
    //
    let abs_err = (value_1.copy() - value_2).norm();
    let rel_err = abs_err / value_2.norm();
    if abs_err < ABS_TOL || rel_err < REL_TOL {
        Ok(())
    } else {
        Err(TestError {
            message: format!(
            "\n\x1b[1;91mAssertion `left ≈= right` failed.\n\x1b[0;91m  left: {}\n right: {}\x1b[0",
            value_1, value_2
        ),
        })
    }
}

pub struct TestError {
    pub message: String,
}

impl fmt::Debug for TestError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\n\x1b[0;2;31m{}\x1b[0m\n",
            self.message,
            get_defeat_message()
        )
    }
}
