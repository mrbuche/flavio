//! Cohesive constitutive models.

#[cfg(test)]
pub mod test;

mod smith_ferrante;

pub use smith_ferrante::SmithFerrante;

use crate::
{
    math::
    {
        TensorRank1Trait,
        TensorRank2Trait
    },
    mechanics::
    {
        Displacement,
        Normal,
        Scalar,
        Stiffness,
        Traction,
        IDENTITY
    }
};
use super::
{
    Constitutive,
    Parameters
};

// TODO: loading/unloading rules using an ISV

/// Required methods for cohesive constitutive models.
pub trait Cohesive<'a>
where
    Self: Constitutive<'a>
{
    fn calculate_traction(&self, displacement: &Displacement, normal: &Normal) -> Traction;
    fn calculate_stiffnesses(&self, displacement: &Displacement, normal: &Normal) -> (Stiffness, Stiffness);
}
