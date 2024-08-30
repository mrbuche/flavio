#[cfg(test)]
mod test;

pub mod list;

use super::Tensor;

/// A tensor of rank 0 (a scalar).
pub type TensorRank0 = f64;

impl Tensor for TensorRank0 {
    fn normy(&self) -> TensorRank0 {
        *self
    }
    fn zeroy() -> Self {
        0.0
    }
}
