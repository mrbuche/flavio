//! Mathematics library.

/// Special mathematical functions.
pub mod special;

mod tensor;

pub const ONE_SIXTH: TensorRank0 = 1.0 / 6.0;
pub const TWO_THIRDS: TensorRank0 = 2.0 / 3.0;
pub const FOUR_THIRDS: TensorRank0 = 4.0 / 3.0;
pub const FIVE_THIRDS: TensorRank0 = 5.0 / 3.0;
pub const SEVEN_THIRDS: TensorRank0 = 7.0 / 3.0;
pub const ONE_TWENTY_FOURTH: TensorRank0 = 1.0 / 24.0;

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
        zero as tensor_rank_1_zero,
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
