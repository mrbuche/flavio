#[cfg(test)]
mod test;

use super::*;

/// A composite hyperelastic constitutive model constructed using the additive decomposition.
pub struct CompositeHyperelasticConstitutiveModelAdditiveDecomposition<C1, C2>
{
    hyperelastic_constitutive_model_1: C1,
    hyperelastic_constitutive_model_2: C2
}

impl<'a, C1, C2> ConstitutiveModel<'a> for CompositeHyperelasticConstitutiveModelAdditiveDecomposition<C1, C2>
where
    C1: ConstitutiveModel<'a>,
    C2: ConstitutiveModel<'a>
{
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma}(\mathbf{F}) = \boldsymbol{\sigma}_1(\mathbf{F}) + \boldsymbol{\sigma}_2(\mathbf{F})
    /// ```
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        self.get_constitutive_model_1().calculate_cauchy_stress(deformation_gradient) + self.get_constitutive_model_2().calculate_cauchy_stress(deformation_gradient)
    }
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\mathcal{T}}(\mathbf{F}) = \boldsymbol{\mathcal{T}}_1(\mathbf{F}) + \boldsymbol{\mathcal{T}}_2(\mathbf{F})
    /// ```
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        self.get_constitutive_model_1().calculate_cauchy_tangent_stiffness(deformation_gradient) + self.get_constitutive_model_2().calculate_cauchy_tangent_stiffness(deformation_gradient)
    }
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a(\mathbf{F}) = a_1(\mathbf{F}) + a_2(\mathbf{F})
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar
    {
        self.get_constitutive_model_1().calculate_helmholtz_free_energy_density(deformation_gradient) + self.get_constitutive_model_2().calculate_helmholtz_free_energy_density(deformation_gradient)
    }
    /// Dummy method that will panic, use [CompositeConstitutiveModel::construct()] instead.
    fn new(_parameters: ConstitutiveModelParameters<'a>) -> Self
    {
        panic!()
    }
}

/// Hyperelastic constitutive model implementation of a composite hyperelastic constitutive model constructed using the additive decomposition.
impl<C1, C2> HyperelasticConstitutiveModel for CompositeHyperelasticConstitutiveModelAdditiveDecomposition<C1, C2>
{
    fn get_bulk_modulus(&self) -> &Scalar
    {
        panic!()
    }
    fn get_shear_modulus(&self) -> &Scalar
    {
        panic!()
    }
}

/// Composite constitutive model implementation of a composite hyperelastic constitutive model constructed using the additive decomposition.
impl<C1, C2> CompositeConstitutiveModel<C1, C2> for CompositeHyperelasticConstitutiveModelAdditiveDecomposition<C1, C2>
{
    fn construct(hyperelastic_constitutive_model_1: C1, hyperelastic_constitutive_model_2: C2) -> Self
    {
        Self
        {
            hyperelastic_constitutive_model_1,
            hyperelastic_constitutive_model_2
        }
    }
    fn get_constitutive_model_1(&self) -> &C1
    {
        &self.hyperelastic_constitutive_model_1
    }
    fn get_constitutive_model_2(&self) -> &C2
    {
        &self.hyperelastic_constitutive_model_2
    }
}
