//! Elastic constitutive models.
//!
//! Elastic constitutive models cannot be defined by a Helmholtz free energy density but still depend on only the deformation gradient.
//! These constitutive models are therefore defined by a relation for some stress measure as a function of the deformation gradient.
//! Consequently, the tangent stiffness associated with the first Piola-Kirchoff stress is not symmetric for elastic constitutive models.
//!
//! ```math
//! \mathcal{C}_{iJkL} \neq \mathcal{C}_{kLiJ}
//! ```

#[cfg(test)]
pub mod test;

mod almansi_hamel;

pub use self::
{
    almansi_hamel::AlmansiHamelModel
};
use super::*;

/// Required methods for elastic constitutive models.
pub trait ElasticConstitutiveModel
{
    /// Returns the bulk modulus.
    fn get_bulk_modulus(&self) -> &Scalar;
    /// Returns the shear modulus.
    fn get_shear_modulus(&self) -> &Scalar;
}
