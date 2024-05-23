//! Mathematics library.

/// Special mathematical functions.
pub mod special;

mod tensor;

pub use tensor::
{
    Convert,
    rank_0::
    {
        TensorRank0,
        list::
        {
            TensorRank0List,
            TensorRank0ListTrait
        }
    },
    rank_1::
    {
        TensorRank1,
        TensorRank1Trait,
        list::
        {
            TensorRank1List,
            TensorRank1ListTrait
        },
        list_2d::
        {
            TensorRank1List2D,
            TensorRank1List2DTrait
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
        },
        list_2d::
        {
            TensorRank2List2D,
            TensorRank2List2DTrait
        }
    },
    rank_3::
    {
        TensorRank3,
        TensorRank3Trait,
        levi_civita,
        list::
        {
            TensorRank3List,
            TensorRank3ListTrait
        },
        list_2d::
        {
            TensorRank3List2D,
            TensorRank3List2DTrait
        },
        list_3d::
        {
            TensorRank3List3D,
            TensorRank3List3DTrait
        }
    },
    rank_4::
    {
        ContractAllIndicesWithFirstIndicesOf,
        ContractFirstSecondIndicesWithSecondIndicesOf,
        ContractFirstThirdFourthIndicesWithFirstIndicesOf,
        ContractSecondIndexWithFirstIndexOf,
        ContractSecondFourthIndicesWithFirstIndicesOf,
        ContractThirdFourthIndicesWithFirstSecondIndicesOf,
        TensorRank4,
        TensorRank4Inverse,
        TensorRank4Trait,
        list::
        {
            TensorRank4List,
            TensorRank4ListTrait
        }
    }
};