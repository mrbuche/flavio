use super::TensorRank0;
use crate::{get_defeat_message, ABS_TOL, EPSILON, REL_TOL};
use std::{cmp::PartialEq, fmt};

#[cfg(test)]
pub trait TensorError {
    fn error(
        &self,
        comparator: &Self,
        tol_abs: &TensorRank0,
        tol_rel: &TensorRank0,
    ) -> Option<usize>;
    fn error_fd(&self, comparator: &Self, epsilon: &TensorRank0) -> Option<usize>;
}

pub fn assert_eq<'a, T: fmt::Display + PartialEq + TensorError>(
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

pub fn assert_eq_from_fd<'a, T: fmt::Display + TensorError>(
    value: &'a T,
    value_fd: &'a T,
) -> Result<(), TestError> {
    if let Some(error_count) = value.error_fd(value_fd, &EPSILON) {
        Err(TestError {
            message: format!(
            "\n\x1b[1;91mAssertion `left ≈= right` failed in {} places.\n\x1b[0;91m  left: {}\n right: {}\x1b[0",
            error_count, value, value_fd
        ),
        })
    } else {
        Ok(())
    }
}

pub fn assert_eq_within_tols<'a, T: fmt::Display + TensorError>(
    value_1: &'a T,
    value_2: &'a T,
) -> Result<(), TestError> {
    if let Some(error_count) = value_1.error(value_2, &ABS_TOL, &REL_TOL) {
        Err(TestError {
            message: format!(
            "\n\x1b[1;91mAssertion `left ≈= right` failed in {} places.\n\x1b[0;91m  left: {}\n right: {}\x1b[0",
            error_count, value_1, value_2
        ),
        })
    } else {
        Ok(())
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
