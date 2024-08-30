use crate::math::{TensorRank0, TensorRank0List};
use std::f64::consts::TAU;

pub fn zero_to_tau<const W: usize>() -> TensorRank0List<W> {
    (0..W)
        .map(|i| TAU * (i as TensorRank0) / ((W - 1) as TensorRank0))
        .collect()
}
