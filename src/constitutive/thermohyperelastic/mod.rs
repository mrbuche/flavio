#[cfg(test)]
mod test;

mod saint_venant_kirchoff;

pub use self::saint_venant_kirchoff::SaintVenantKirchoffModel;
use super::
{
    *, thermoelastic::ThermoelasticConstitutiveModel
};

/// Required methods for thermohyperelastic constitutive models.
pub trait ThermohyperelasticConstitutiveModel
where
    Self: ThermoelasticConstitutiveModel
{
    /// Calculates and returns the Helmholtz free energy density.
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient, temperature: &Temperature) -> Scalar;
}
