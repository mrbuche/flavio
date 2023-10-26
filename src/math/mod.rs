mod tensor;

pub use tensor::
{
    rank_0::TensorRank0,
    rank_1::
    {
        TensorRank1,
        TensorRank1Trait,
        list::
        {
            TensorRank1List,
            TensorRank1ListTrait
        }
    },
    rank_2::
    {
        TensorRank2,
        TensorRank2Trait,
        list::
        {
            TensorRank2List,
            TensorRank2ListTrait
        }
    },
    rank_3::
    {
        TensorRank3,
        TensorRank3Trait
    },
    rank_4::
    {
        TensorRank4,
        TensorRank4Trait
    }
};