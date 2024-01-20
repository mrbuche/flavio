//! Elastic constitutive models.
//!
//! Elastic constitutive models cannot be derived from a Helmholtz free energy density, they are instead defined by a relation for the stress. As a result, the tangent stiffness associated with the first Piola-Kirchoff stress is not symmetric for elastic constitutive models.
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
