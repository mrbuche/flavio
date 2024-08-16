#[cfg(test)]
mod test;

use crate::
{
    constitutive::
    {
        ConstitutiveError,
        hybrid::
        {
            Additive,
            Hybrid
        },
        solid::hyperelastic::Hyperelastic
    },
    mechanics::
    {
        DeformationGradient,
        Scalar
    }
};

/// Hyperelastic constitutive model implementation of hybrid hyperelastic constitutive models that are based on the additive decomposition.
impl<'a, C1: Hyperelastic<'a>, C2: Hyperelastic<'a>> Hyperelastic<'a> for Additive<C1, C2>
{
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a(\mathbf{F}) = a_1(\mathbf{F}) + a_2(\mathbf{F})
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Result<Scalar, ConstitutiveError>
    {
        Ok(self.get_constitutive_model_1().calculate_helmholtz_free_energy_density(deformation_gradient)? + self.get_constitutive_model_2().calculate_helmholtz_free_energy_density(deformation_gradient)?)
    }
}
