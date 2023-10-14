#[cfg(test)]
mod test;

use super::
{
    rank_0::TensorRank0,
    // rank_2::TensorRank2
};

#[derive(Clone, Copy)]
pub struct TensorRank1<const D: usize>
(
    pub [TensorRank0; D]
);
