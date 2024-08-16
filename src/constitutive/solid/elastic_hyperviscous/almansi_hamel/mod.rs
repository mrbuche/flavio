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
pub struct AlmansiHamel<'a> {
    parameters: Parameters<'a>,
}

/// Constitutive model implementation of the Almansi-Hamel viscoelastic constitutive model.
impl<'a> Constitutive<'a> for AlmansiHamel<'a> {
    fn new(parameters: Parameters<'a>) -> Self {
        Self { parameters }
    }
}

/// Solid constitutive model implementation of the Almansi-Hamel viscoelastic constitutive model.
impl<'a> Solid<'a> for AlmansiHamel<'a> {
    fn get_bulk_modulus(&self) -> &Scalar {
        &self.parameters[0]
    }
    fn get_shear_modulus(&self) -> &Scalar {
        &self.parameters[1]
    }
}

/// Viscous constitutive model implementation of the Almansi-Hamel viscoelastic constitutive model.
impl<'a> Viscous<'a> for AlmansiHamel<'a> {
    fn get_bulk_viscosity(&self) -> &Scalar {
        &self.parameters[2]
    }
    fn get_shear_viscosity(&self) -> &Scalar {
        &self.parameters[3]
    }
}

/// Viscoelastic constitutive model implementation of the Almansi-Hamel viscoelastic constitutive model.
impl<'a> Viscoelastic<'a> for AlmansiHamel<'a> {
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \mathbf{}(\mathbf{F},\dot\mathbf{F}) = 2\mu\mathbf{e}' + \kappa\,\mathrm{tr}(\mathbf{e})\mathbf{1} + 2\eta\mathbf{D}' + \zeta\,\mathrm{tr}(\mathbf{D})\mathbf{1}
    /// ```
    fn calculate_cauchy_stress(
        &self,
        deformation_gradient: &DeformationGradient,
        deformation_gradient_rate: &DeformationGradientRate,
    ) -> CauchyStress {
        let (inverse_deformation_gradient, jacobian) =
            deformation_gradient.inverse_and_determinant();
        let strain = (IDENTITY
            - inverse_deformation_gradient.transpose() * &inverse_deformation_gradient)
            * 0.5;
        let (deviatoric_strain, strain_trace) = strain.deviatoric_and_trace();
        let velocity_gradient = deformation_gradient_rate * inverse_deformation_gradient;
        let strain_rate = (&velocity_gradient + velocity_gradient.transpose()) * 0.5;
        let (deviatoric_strain_rate, strain_rate_trace) = strain_rate.deviatoric_and_trace();
        deviatoric_strain * (2.0 * self.get_shear_modulus() / jacobian)
            + deviatoric_strain_rate * (2.0 * self.get_shear_viscosity() / jacobian)
            + IDENTITY
                * ((self.get_bulk_modulus() * strain_trace
                    + self.get_bulk_viscosity() * strain_rate_trace)
                    / jacobian)
    }
    /// Calculates and returns the rate tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{V}_{IJkL}(\mathbf{F}) = \eta\,\delta_{ik}F_{jL}^{-T} + \eta\,\delta_{jk}F_{iL}^{-T} + \left(\zeta - \frac{2}{3}\,\eta\right)\delta_{ij}F_{kL}^{-T}
    /// ```
    fn calculate_cauchy_rate_tangent_stiffness(
        &self,
        deformation_gradient: &DeformationGradient,
        _: &DeformationGradientRate,
    ) -> CauchyRateTangentStiffness {
        let (deformation_gradient_inverse_transpose, jacobian) =
            deformation_gradient.inverse_transpose_and_determinant();
        let scaled_deformation_gradient_inverse_transpose =
            &deformation_gradient_inverse_transpose * self.get_shear_viscosity() / jacobian;
        CauchyRateTangentStiffness::dyad_ik_jl(
            &IDENTITY,
            &scaled_deformation_gradient_inverse_transpose,
        ) + CauchyRateTangentStiffness::dyad_il_jk(
            &scaled_deformation_gradient_inverse_transpose,
            &IDENTITY,
        ) + CauchyRateTangentStiffness::dyad_ij_kl(
            &(IDENTITY
                * ((self.get_bulk_viscosity() - TWO_THIRDS * self.get_shear_viscosity())
                    / jacobian)),
            &deformation_gradient_inverse_transpose,
        )
    }
}

/// Elastic-hyperviscous constitutive model implementation of the Almansi-Hamel elastic-hyperviscous constitutive model.
impl<'a> ElasticHyperviscous<'a> for AlmansiHamel<'a> {
    /// Calculates and returns the viscous dissipation.
    ///
    /// ```math
    /// \phi(\mathbf{F},\dot{\mathbf{F}}) = \eta\,\mathrm{tr}(\mathbf{D}^2) + \frac{1}{2}\left(\zeta - \frac{2}{3}\,\eta\right)\mathrm{tr}(\mathbf{D})^2
    /// ```
    fn calculate_viscous_dissipation(
        &self,
        deformation_gradient: &DeformationGradient,
        deformation_gradient_rate: &DeformationGradientRate,
    ) -> Scalar {
        let velocity_gradient = deformation_gradient_rate * deformation_gradient.inverse();
        let strain_rate = (&velocity_gradient + velocity_gradient.transpose()) * 0.5;
        self.get_shear_viscosity() * strain_rate.squared_trace()
            + 0.5
                * (self.get_bulk_viscosity() - TWO_THIRDS * self.get_shear_viscosity())
                * strain_rate.trace().powi(2)
    }
}
