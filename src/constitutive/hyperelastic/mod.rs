#[cfg(test)]
mod test;

mod additive_decomposition;
mod arruda_boyce;
mod gent;
mod mooney_rivlin;
mod neo_hookean;
mod yeoh;

pub use self::
{
    additive_decomposition::CompositeHyperelasticConstitutiveModelAdditiveDecomposition,
    arruda_boyce::ArrudaBoyceModel,
    gent::GentModel,
    mooney_rivlin::MooneyRivlinModel,
    neo_hookean::NeoHookeanModel,
    yeoh::YeohModel
};
use super::*;

/// Required methods for hyperelastic constitutive models.
pub trait HyperelasticConstitutiveModel
{
    /// Returns the bulk modulus.
    fn get_bulk_modulus(&self) -> &Scalar;
    /// Returns the shear modulus.
    fn get_shear_modulus(&self) -> &Scalar;
}
