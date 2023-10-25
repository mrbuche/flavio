mod tensor;

pub use tensor::
{
    rank_0::TensorRank0,
    rank_1::
    {
        TensorRank1,
        TensorRank1Traits,
        list::
        {
            TensorRank1List,
            TensorRank1ListTraits
        }
    },
    rank_2::
    {
        TensorRank2,
        TensorRank2Traits
    },
    rank_3::
    {
        TensorRank3,
        TensorRank3Traits
    },
    rank_4::
    {
        TensorRank4,
        TensorRank4Traits
    }
};