//! Elastoplastic constitutive models.

#[cfg(test)]
pub mod test;

use super::
{
    elastic::Elastic,
    plastic::Plastic
};

/// Required methods for elastoplastic constitutive models.
pub trait Elastoplastic<'a>
where
    Self: Elastic<'a> + Plastic<'a>
{
    // maybe once Fp is an ISV you can make a "calculate_Fe" method
}
