//! Hyperelastic constitutive models.
//!
//! Hyperelastic constitutive models are completely defined by a Helmholtz free energy density function of the deformation gradient.
//! Consequently, the tangent stiffness associated with the first Piola-Kirchoff stress is symmetric for hyperelastic constitutive models.
//!
//! ```math
//! \mathcal{C}_{iJkL} = \mathcal{C}_{kLiJ}
//! ```

#[cfg(test)]
pub mod test;

mod additive_decomposition;
mod arruda_boyce;
mod fung;
mod gent;
mod mooney_rivlin;
mod neo_hookean;
mod saint_venant_kirchoff;
mod yeoh;

pub use self::
{
    arruda_boyce::ArrudaBoyce,
    fung::Fung,
    gent::Gent,
    mooney_rivlin::MooneyRivlin,
    neo_hookean::NeoHookean,
    saint_venant_kirchoff::SaintVenantKirchoff,
    yeoh::Yeoh
};
use super::
{
    *, elastic::ElasticConstitutiveModel
};

/// Required methods for hyperelastic constitutive models.
pub trait HyperelasticConstitutiveModel<'a>
where
    Self: ElasticConstitutiveModel<'a>
{
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a = a(\mathbf{F})
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar;
}

/// A composite hyperelastic constitutive model.
pub struct CompositeHyperelasticConstitutiveModel<C1, C2>
{
    hyperelastic_constitutive_model_1: C1,
    hyperelastic_constitutive_model_2: C2
}
