//! Fluid constitutive models.

pub mod plastic;
pub mod viscous;

use crate::
{
    math::TensorRank2Trait,
    mechanics::
    {
        MandelStress,
        Scalar,
        StretchingRatePlastic
    }
};
use super::Constitutive;

/// Required methods for fluid constitutive models.
pub trait Fluid<'a>
where
    Self: Constitutive<'a>
{}
