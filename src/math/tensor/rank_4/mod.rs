#[cfg(test)]
mod test;

use super::
{
    rank_2::TensorRank2
};

pub struct TensorRank4<const D: usize>
(
    pub [[TensorRank2<D>; D]; D]
);
