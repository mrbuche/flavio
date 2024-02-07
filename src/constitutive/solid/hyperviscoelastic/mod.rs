//! Hyperviscoelastic constitutive models.

#[cfg(test)]
mod test;

mod saint_venant_kirchoff;

pub use saint_venant_kirchoff::SaintVenantKirchoff;

use super::
{
    *,
    viscoelastic::Viscoelastic,
    super::fluid::viscous::Viscous
};

/// Required methods for hyperviscoelastic constitutive models.
pub trait Hyperviscoelastic<'a>
where
    Self: Viscoelastic<'a>
{
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a = a(\mathbf{F})
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar;
    /// Calculates and returns the viscous dissipation.
    ///
    /// ```math
    /// b = b(\dot{\mathbf{F}})
    /// ```
    fn calculate_viscous_dissipation(&self, deformation_gradient_dot: &DeformationGradientDot) -> Scalar;
}
