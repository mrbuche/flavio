//! Elastic-hyperviscous constitutive models.
//!
//! Elastic-hyperviscous constitutive models are defined by an elastic stress tensor function and a viscous dissipation function.
//!
//! ```math
//! \mathbf{P}:\dot{\mathbf{F}} - \mathbf{P}^e(\mathbf{F}):\dot{\mathbf{F}} - \phi(\mathbf{F},\dot{\mathbf{F}}) \geq 0
//! ```
//! Satisfying the second law of thermodynamics though a minimum viscous dissipation principal yields a relation for the stress.
//!
//! ```math
//! \mathbf{P} = \mathbf{P}^e + \frac{\partial\phi}{\partial\dot{\mathbf{F}}}
//! ```
//! Consequently, the rate tangent stiffness associated with the first Piola-Kirchoff stress is symmetric for elastic-hyperviscous models.
//!
//! ```math
//! \mathcal{U}_{iJkL} = \mathcal{U}_{kLiJ}
//! ```

#[cfg(test)]
pub mod test;

mod almansi_hamel;

pub use almansi_hamel::AlmansiHamel;

use super::
{
    *,
    viscoelastic::Viscoelastic,
    super::fluid::viscous::Viscous
};

/// Required methods for elastic-hyperviscous constitutive models.
pub trait ElasticHyperviscous<'a>
where
    Self: Viscoelastic<'a>
{
    /// Calculates and returns the dissipation potential.
    ///
    /// ```math
    /// \mathbf{P}^e(\mathbf{F}):\dot{\mathbf{F}} + \phi(\mathbf{F},\dot{\mathbf{F}})
    /// ```
    fn calculate_dissipation_potential(&self, deformation_gradient: &DeformationGradient, deformation_gradient_rate: &DeformationGradientRate) -> Scalar
    {
        self.calculate_first_piola_kirchoff_stress(deformation_gradient, &DeformationGradientRate::zero()).full_contraction(deformation_gradient_rate) + self.calculate_viscous_dissipation(deformation_gradient, deformation_gradient_rate)
    }
    /// Calculates and returns the viscous dissipation.
    ///
    /// ```math
    /// \phi = \phi(\mathbf{F},\dot{\mathbf{F}})
    /// ```
    fn calculate_viscous_dissipation(&self, deformation_gradient: &DeformationGradient, deformation_gradient_rate: &DeformationGradientRate) -> Scalar;
}
