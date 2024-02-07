#[cfg(test)]
mod test;

use super::*;

impl<'a, C1, C2> Constitutive<'a> for CombinedHyperelastic<C1, C2>
where
    C1: Hyperelastic<'a>,
    C2: Hyperelastic<'a>
{
    /// Dummy method that will panic, use [Self::construct()] instead.
    fn new(_parameters: Parameters<'a>) -> Self
    {
        panic!()
    }
}

impl<'a, C1, C2> Solid<'a> for CombinedHyperelastic<C1, C2>
where
    C1: Hyperelastic<'a>,
    C2: Hyperelastic<'a>
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

impl<'a, C1, C2> Elastic<'a> for CombinedHyperelastic<C1, C2>
where
    C1: Hyperelastic<'a>,
    C2: Hyperelastic<'a>
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
}

/// Hyperelastic constitutive model implementation of a combined hyperelastic constitutive model constructed using the additive decomposition.
impl<'a, C1, C2> Hyperelastic<'a> for CombinedHyperelastic<C1, C2>
where
    C1: Hyperelastic<'a>,
    C2: Hyperelastic<'a>
{
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a(\mathbf{F}) = a_1(\mathbf{F}) + a_2(\mathbf{F})
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar
    {
        self.get_constitutive_model_1().calculate_helmholtz_free_energy_density(deformation_gradient) + self.get_constitutive_model_2().calculate_helmholtz_free_energy_density(deformation_gradient)
    }
}

/// Combined constitutive model implementation of a combined hyperelastic constitutive model constructed using the additive decomposition.
impl<'a, C1, C2> Combined<C1, C2> for CombinedHyperelastic<C1, C2>
where
    C1: Hyperelastic<'a>,
    C2: Hyperelastic<'a>
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
