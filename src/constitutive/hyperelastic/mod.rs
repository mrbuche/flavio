//! Hyperelastic constitutive models.
//!
//! Hyperelastic constitutive models are completely defined by a Helmholtz free energy density function of the deformation gradient.
//! Consequently, the tangent stiffness associated with the first Piola-Kirchoff stress is symmetric for hyperelastic constitutive models.
//!
//! ```math
//! \mathcal{C}_{iJkL} = \mathcal{C}_{kLiJ}
//! ```

#[cfg(test)]
mod test;

// mod additive_decomposition;
// mod arruda_boyce;
// mod fung;
// mod gent;
// mod mooney_rivlin;
mod neo_hookean;
// mod saint_venant_kirchoff;
// mod yeoh;

pub use self::
{
    // arruda_boyce::ArrudaBoyceModel,
    // fung::FungModel,
    // gent::GentModel,
    // mooney_rivlin::MooneyRivlinModel,
    neo_hookean::NeoHookeanModel,
    // saint_venant_kirchoff::SaintVenantKirchoffModel,
    // yeoh::YeohModel
};
use super::
{
    *, elastic::ElasticConstitutiveModel
};

/// Required methods for hyperelastic constitutive models.
pub trait HyperelasticConstitutiveModel
where
    Self: ElasticConstitutiveModel
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
