//! Hyperelastic constitutive models.
//!
//! Hyperelastic constitutive models are derived from and defined by a helmholtz free energy density. As a result, the tangent stiffness associated with the first Piola-Kirchoff stress is symmetric for hyperelastic constitutive models.
//!
//! ```math
//! \mathcal{C}_{iJkL} = \mathcal{C}_{kLiJ}
//! ```

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
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a = a(\mathbf{F})
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar;
    /// Returns the bulk modulus.
    fn get_bulk_modulus(&self) -> &Scalar;
    /// Returns the shear modulus.
    fn get_shear_modulus(&self) -> &Scalar;
}

/// A composite hyperelastic constitutive model.
pub struct CompositeHyperelasticConstitutiveModel<C1, C2>
{
    hyperelastic_constitutive_model_1: C1,
    hyperelastic_constitutive_model_2: C2
}