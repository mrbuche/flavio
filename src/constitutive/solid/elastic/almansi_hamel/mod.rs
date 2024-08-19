#[cfg(test)]
mod test;

use super::*;

/// The Almansi-Hamel elastic constitutive model.
///
/// **Parameters**
/// - The bulk modulus $`\kappa`$.
/// - The shear modulus $`\mu`$.
///
/// **External variables**
/// - The deformation gradient $`\mathbf{F}`$.
///
/// **Internal variables**
/// - None.
///
/// **Notes**
/// - The Almansi-Hamel strain measure is given by $`\mathbf{e}=\tfrac{1}{2}(\mathbf{1}-\mathbf{B}^{-1})`$.
#[derive(Debug)]
pub struct AlmansiHamel<'a> {
    parameters: Parameters<'a>,
}

/// Constitutive model implementation of the Almansi-Hamel elastic constitutive model.
impl<'a> Constitutive<'a> for AlmansiHamel<'a> {
    fn new(parameters: Parameters<'a>) -> Self {
        Self { parameters }
    }
}

/// Solid constitutive model implementation of the Almansi-Hamel elastic constitutive model.
impl<'a> Solid<'a> for AlmansiHamel<'a> {
    fn get_bulk_modulus(&self) -> &Scalar {
        &self.parameters[0]
    }
    fn get_shear_modulus(&self) -> &Scalar {
        &self.parameters[1]
    }
}

/// Elastic constitutive model implementation of the Almansi-Hamel elastic constitutive model.
impl<'a> Elastic<'a> for AlmansiHamel<'a> {
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma}(\mathbf{F}) = \frac{2\mu}{J}\,\mathbf{e}' + \frac{\kappa}{J}\,\mathrm{tr}(\mathbf{e})\mathbf{1}
    /// ```
    fn calculate_cauchy_stress(
        &self,
        deformation_gradient: &DeformationGradient,
    ) -> Result<CauchyStress, ConstitutiveError> {
        let jacobian = deformation_gradient.determinant();
        if jacobian > 0.0 {
            let inverse_deformation_gradient = deformation_gradient.inverse();
            let strain = (IDENTITY
                - inverse_deformation_gradient.transpose() * &inverse_deformation_gradient)
                * 0.5;
            let (deviatoric_strain, strain_trace) = strain.deviatoric_and_trace();
            Ok(
                deviatoric_strain * (2.0 * self.get_shear_modulus() / jacobian)
                    + IDENTITY * (self.get_bulk_modulus() * strain_trace / jacobian),
            )
        } else {
            Err(ConstitutiveError::InvalidJacobian(
                jacobian,
                deformation_gradient.copy(),
                format!("{:?}", &self),
            ))
        }
    }
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{T}_{ijkL}(\mathbf{F}) = \frac{\mu}{J}\left[B_{jk}^{-1}F_{iL}^{-T} + B_{ik}^{-1}F_{jL}^{-T} - \frac{2}{3}\,\delta_{ij}B_{km}^{-1}F_{mL}^{-T} - 2e_{ij}'F_{kL}^{-T}\right] + \frac{\kappa}{J}\left[\delta_{ij}B_{km}^{-1}F_{mL}^{-T} - \mathrm{tr}(\mathbf{e})\delta_{ij}F_{kL}^{-T}\right]
    /// ```
    fn calculate_cauchy_tangent_stiffness(
        &self,
        deformation_gradient: &DeformationGradient,
    ) -> Result<CauchyTangentStiffness, ConstitutiveError> {
        let jacobian = deformation_gradient.determinant();
        if jacobian > 0.0 {
            let inverse_transpose_deformation_gradient = deformation_gradient.inverse_transpose();
            let inverse_left_cauchy_green_deformation = &inverse_transpose_deformation_gradient
                * inverse_transpose_deformation_gradient.transpose();
            let strain = (IDENTITY - &inverse_left_cauchy_green_deformation) * 0.5;
            let (deviatoric_strain, strain_trace) = strain.deviatoric_and_trace();
            Ok((CauchyTangentStiffness::dyad_il_jk(
                &inverse_transpose_deformation_gradient,
                &inverse_left_cauchy_green_deformation,
            ) + CauchyTangentStiffness::dyad_ik_jl(
                &inverse_left_cauchy_green_deformation,
                &inverse_transpose_deformation_gradient,
            )) * (self.get_shear_modulus() / jacobian)
                + CauchyTangentStiffness::dyad_ij_kl(
                    &IDENTITY,
                    &(inverse_left_cauchy_green_deformation
                        * &inverse_transpose_deformation_gradient
                        * ((self.get_bulk_modulus() - self.get_shear_modulus() * TWO_THIRDS)
                            / jacobian)),
                )
                - CauchyTangentStiffness::dyad_ij_kl(
                    &(deviatoric_strain * (2.0 * self.get_shear_modulus() / jacobian)
                        + IDENTITY * (self.get_bulk_modulus() * strain_trace / jacobian)),
                    &inverse_transpose_deformation_gradient,
                ))
        } else {
            Err(ConstitutiveError::InvalidJacobian(
                jacobian,
                deformation_gradient.copy(),
                format!("{:?}", &self),
            ))
        }
    }
}
