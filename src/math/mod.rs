//! Mathematics library.

/// Special mathematical functions.
pub mod special;

mod tensor;

pub const FOUR_THIRDS: TensorRank0 = 4.0 / 3.0;
pub const FIVE_THIRDS: TensorRank0 = 5.0 / 3.0;
pub const NINE_FIFTHS: TensorRank0 = 9.0 / 5.0;
pub const ONE_SIXTH: TensorRank0 = 1.0 / 6.0;
pub const ONE_TWENTY_FOURTH: TensorRank0 = 1.0 / 24.0;
pub const SEVEN_THIRDS: TensorRank0 = 7.0 / 3.0;
pub const TWO_THIRDS: TensorRank0 = 2.0 / 3.0;

pub use tensor::{
    rank_0::{
        list::{TensorRank0List, TensorRank0ListTrait},
        TensorRank0,
    },
    rank_1::{
        list::{TensorRank1List, TensorRank1ListTrait},
        list_2d::{TensorRank1List2D, TensorRank1List2DTrait},
        zero as tensor_rank_1_zero, TensorRank1, TensorRank1Trait,
    },
    rank_2::{
        list::{TensorRank2List, TensorRank2ListTrait},
        list_2d::{TensorRank2List2D, TensorRank2List2DTrait},
        TensorRank2, TensorRank2Trait,
    },
    rank_3::{
        levi_civita,
        list::{TensorRank3List, TensorRank3ListTrait},
        list_2d::{TensorRank3List2D, TensorRank3List2DTrait},
        list_3d::{TensorRank3List3D, TensorRank3List3DTrait},
        TensorRank3, TensorRank3Trait,
    },
    rank_4::{
        list::{TensorRank4List, TensorRank4ListTrait},
        ContractAllIndicesWithFirstIndicesOf, ContractFirstSecondIndicesWithSecondIndicesOf,
        ContractFirstThirdFourthIndicesWithFirstIndicesOf,
        ContractSecondFourthIndicesWithFirstIndicesOf, ContractSecondIndexWithFirstIndexOf,
        ContractThirdFourthIndicesWithFirstSecondIndicesOf, TensorRank4, TensorRank4Trait,
    },
    Convert,
};
