#[cfg(test)]
mod test;

mod saint_venant_kirchoff;

pub use self::saint_venant_kirchoff::SaintVenantKirchoffModel;
use super::
{
    *,
    hyperelastic::HyperelasticConstitutiveModel
};

/// Required methods for thermohyperelastic constitutive models.
pub trait ThermohyperelasticConstitutiveModel<V>
where
    Self: HyperelasticConstitutiveModel<V>
{
    /// Returns the coefficient of thermal expansion.
    fn get_coefficient_of_thermal_expansion(&self) -> &Scalar;
    /// Returns the reference temperature.
    fn get_reference_temperature(&self) -> &Scalar;
}
