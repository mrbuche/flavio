#[cfg(test)]
mod test;

use super::
{
    rank_1::TensorRank1
};

pub struct TensorRank2<const D: usize>
(
    pub [TensorRank1<D>; D]
);
