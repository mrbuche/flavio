#[cfg(test)]
mod test;

use super::*;

/// The Almansi-Hamel viscoelastic constitutive model.
///
/// **Parameters**
/// - The bulk modulus $`\kappa`$.
/// - The shear modulus $`\mu`$.
/// - The bulk viscosity $`\zeta`$.
/// - The shear viscosity $`\eta`$.
///
/// **External variables**
/// - The deformation gradient $`\mathbf{F}`$.
/// - The deformation gradient rate $`\dot{\mathbf{F}}`$.
///
/// **Internal variables**
/// - None.
///
/// **Notes**
/// - The Almansi-Hamel strain measure is given by $`\mathbf{e}=\tfrac{1}{2}(\mathbf{1}-\mathbf{B}^{-1})`$.
pub struct AlmansiHamel<'a>
{
    parameters: Parameters<'a>
}

/// Constitutive model implementation of the Almansi-Hamel viscoelastic constitutive model.
impl<'a> Constitutive<'a> for AlmansiHamel<'a>
{
    fn new(parameters: Parameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

/// Solid constitutive model implementation of the Almansi-Hamel viscoelastic constitutive model.
impl<'a> Solid<'a> for AlmansiHamel<'a>
{
    fn get_bulk_modulus(&self) -> &Scalar
    {
        &self.parameters[0]
    }
    fn get_shear_modulus(&self) -> &Scalar
    {
        &self.parameters[1]
    }
}

/// Viscous constitutive model implementation of the Almansi-Hamel viscoelastic constitutive model.
impl<'a> Viscous<'a> for AlmansiHamel<'a>
{
    fn get_bulk_viscosity(&self) -> &Scalar
    {
        &self.parameters[2]
    }
    fn get_shear_viscosity(&self) -> &Scalar
    {
        &self.parameters[3]
    }
}

/// Viscoelastic constitutive model implementation of the Almansi-Hamel viscoelastic constitutive model.
impl<'a> Viscoelastic<'a> for AlmansiHamel<'a>
{
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma}(\mathbf{F},\dot{\mathbf{F}}) = \frac{2\mu}{J}\,\mathbf{e}' + \frac{\kappa}{J}\,\mathrm{tr}(\mathbf{e})\mathbf{1} + \frac{2\eta}{J}\,\dot{\mathbf{e}}' + \frac{\zeta}{J}\,\mathrm{tr}(\dot{\mathbf{e}})\mathbf{1}
    /// ```
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient, deformation_gradient_rate: &DeformationGradientRate) -> CauchyStress
    {
        todo!()
    }
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{T}_{ijkL}(\mathbf{F},\dot{\mathbf{F}}) = ?
    /// ```
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient, deformation_gradient_rate: &DeformationGradientRate) -> CauchyTangentStiffness
    {
        todo!()
    }
    /// Calculates and returns the rate tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{V}_{ijkL}(\mathbf{F},\dot{\mathbf{F}}) = ?
    /// ```
    fn calculate_cauchy_rate_tangent_stiffness(&self, deformation_gradient: &DeformationGradient, deformation_gradient_rate: &DeformationGradientRate) -> CauchyTangentStiffness
    {
        todo!()
    }
}