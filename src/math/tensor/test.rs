use super::Tensor;
use crate::{get_defeat_message, ABS_TOL, REL_TOL};
use std::{convert::From, fmt};

pub fn assert_eq_within_tols_new<'a, T: Tensor>(
    value_1: &'a T,
    value_2: &'a T,
) -> Result<(), NotEqualWithinTolsValues<'a, T>> {
    let abs_err = (value_1.copy() - value_2).norm();
    let satisfies_rel_tol = false; // or use abs_err (value_1 / value_2 - T::identity()).norm() < REL_TOL;
    if abs_err < ABS_TOL || satisfies_rel_tol {
        Ok(())
    } else {
        Err(NotEqualWithinTolsValues(&value_1, &value_2))
    }
}

pub struct NotEqualWithinTolsValues<'a, T: Tensor>(&'a T, &'a T);

pub struct NotEqualWithinTols(String);

// Display impl might be better to require for Tensor to implement?

impl<'a, T: Tensor> From<NotEqualWithinTolsValues<'a, T>> for NotEqualWithinTols {
    fn from(error: NotEqualWithinTolsValues<'a, T>) -> NotEqualWithinTols {
        NotEqualWithinTols(format!(
            "\n\x1b[1;91mAssertion `left ≈ right` failed.\n\x1b[0;91m  left: {:?}\n right: {:?}\n\x1b[0;2;31m{}\x1b[0m\n",
            error.0,
            error.1,
            get_defeat_message()
        ))
    }
}

impl fmt::Debug for NotEqualWithinTols {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}
