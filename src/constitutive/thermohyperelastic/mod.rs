#[cfg(test)]
mod test;

mod saint_venant_kirchoff;

pub use self::saint_venant_kirchoff::SaintVenantKirchoffModel;
use super::*;

/// Required methods for thermohyperelastic constitutive models.
pub trait ThermohyperelasticConstitutiveModel
{
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a = a(\mathbf{F}, T)
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient, temperature: &Temperature) -> Scalar;
    /// Returns the bulk modulus.
    fn get_bulk_modulus(&self) -> &Scalar;
    /// Returns the shear modulus.
    fn get_shear_modulus(&self) -> &Scalar;
    /// Returns the coefficient of thermal expansion.
    fn get_coefficient_of_thermal_expansion(&self) -> &Scalar;
    /// Returns the reference temperature.
    fn get_reference_temperature(&self) -> &Scalar;
}
