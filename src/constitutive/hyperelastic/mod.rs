#[cfg(test)]
mod test;

pub mod additive_decomposition;

pub mod arruda_boyce;
pub mod gent;
pub mod mooney_rivlin;
pub mod neo_hookean;
pub mod yeoh;

use super::*;

pub trait HyperelasticConstitutiveModel
{
    fn get_bulk_modulus(&self) -> &Scalar;
    fn get_shear_modulus(&self) -> &Scalar;
}
