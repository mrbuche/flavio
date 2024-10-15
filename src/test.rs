use super::{get_defeat_message, ABS_TOL, REL_TOL};
use std::fmt;

pub fn assert_eq_within_tols(value_1: &f64, value_2: &f64) {
    assert!(check_eq_within_tols(value_1, value_2))
}

pub fn check_eq_within_tols(value_1: &f64, value_2: &f64) -> bool {
    (value_1 - value_2).abs() < ABS_TOL || (value_1 / value_2 - 1.0).abs() < REL_TOL
}

pub fn assert_eq_within_tols_new(value_1: f64, value_2: f64) -> Result<(), Foo> {
    if check_eq_within_tols(&value_1, &value_2) {
        Ok(())
    } else {
        Err(Foo(value_1, value_2))
    }
}

pub struct Foo(f64, f64);

impl fmt::Debug for Foo {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "\n\x1b[1;91mfoo\n\x1b[0;91m  left: {}\n right: {}\n\x1b[0;2;31m{}\x1b[0m\n",
            self.0,
            self.1,
            get_defeat_message()
        )
    }
}
