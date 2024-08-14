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

/// Required methods for hyperviscoelastic constitutive models.
pub trait Hyperviscoelastic<'a>
where
    Self: ElasticHyperviscous<'a>
{
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a = a(\mathbf{F})
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Result<Scalar, MyError>;
}



use crate::get_error_message;

#[derive(Debug)]
pub enum MyError {
    InvalidJacobian(Scalar),
}
// Shouldn't determinant() have Err(Invalid) and pass that up here?

impl std::fmt::Display for MyError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let error = match self {
            Self::InvalidJacobian(jacobian) => format!("Invalid Jacobian of {}", jacobian)
        };
        write!(f, "\x1b[91m{}\n\x1b[0;2;31m{}\x1b[0m\n", error, get_error_message())
    }
}
