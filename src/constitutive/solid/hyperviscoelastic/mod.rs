//! Hyperviscoelastic constitutive models.
//!
//! ---
//!
//! Hyperviscoelastic constitutive models are defined by a Helmholtz free energy density function and a viscous dissipation function.
//!
//! ```math
//! \mathbf{P}:\dot{\mathbf{F}} - \dot{a}(\mathbf{F}) - \phi(\mathbf{F},\dot{\mathbf{F}}) \geq 0
//! ```
//! Satisfying the second law of thermodynamics though a minimum viscous dissipation principal yields a relation for the stress.
//!
//! ```math
//! \mathbf{P} = \frac{\partial a}{\partial\mathbf{F}} + \frac{\partial\phi}{\partial\dot{\mathbf{F}}}
//! ```
//! Consequently, the rate tangent stiffness associated with the first Piola-Kirchoff stress is symmetric for hyperviscoelastic models.
//!
//! ```math
//! \mathcal{U}_{iJkL} = \mathcal{U}_{kLiJ}
//! ```

#[cfg(test)]
pub mod test;

mod saint_venant_kirchoff;

pub use saint_venant_kirchoff::SaintVenantKirchoff;

use super::
{
    *,
    elastic_hyperviscous::ElasticHyperviscous,
    viscoelastic::Viscoelastic,
    super::fluid::viscous::Viscous
};
use std::fmt::Debug;

/// Required methods for hyperviscoelastic constitutive models.
pub trait Hyperviscoelastic<'a>
where
    Self: Debug + ElasticHyperviscous<'a>
{
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a = a(\mathbf{F})
    /// ```
    fn calculate_helmholtz_free_energy_density(&'a self, deformation_gradient: &'a DeformationGradient) -> Result<Scalar, ConstitutiveError>
    {
        let jacobian = deformation_gradient.determinant();
        if jacobian > 0.0 {
            Ok(self.calculate_helmholtz_free_energy_density_inner(deformation_gradient))
        } else {
            Err(ConstitutiveError::InvalidJacobian(jacobian, deformation_gradient, format!("{:?}", &self)))
        }
    }
    #[doc(hidden)]
    fn calculate_helmholtz_free_energy_density_inner(&self, deformation_gradient: &DeformationGradient) -> Scalar;

    this is good but you lost the docstring in SVK!

}
