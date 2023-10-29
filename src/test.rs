use super::
{
    ABS_TOL,
    REL_TOL
};

pub fn assert_eq_within_tols(value_1: &f64, value_2: &f64)
{
    assert!(
        (value_1 - value_2).abs() < ABS_TOL ||
        (value_1/value_2 - 1.0).abs() < REL_TOL
    )
}
