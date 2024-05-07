#[cfg(test)]
pub mod test;

mod smith_ferrante;

use crate::
{
    math::TensorRank1Trait,
    mechanics::
    {
        Displacement,
        Normal,
        Scalar,
        Traction
    }
};
use super::
{
    Constitutive,
    Parameters
};

// TODO: loading/unloading rules using an ISV

pub trait Cohesive<'a>
where
    Self: Constitutive<'a>
{
    fn calculate_traction(&self, displacement: &Displacement, normal: &Normal) -> Traction;
}