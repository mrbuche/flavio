//! Elastic constitutive models.
//!
//! ---
//!
//! Elastic constitutive models cannot be defined by a Helmholtz free energy density but still depend on only the deformation gradient.
//! These constitutive models are therefore defined by a relation for some stress measure as a function of the deformation gradient.
//! Consequently, the tangent stiffness associated with the first Piola-Kirchoff stress is not symmetric for elastic constitutive models.
//!
//! ```math
//! \mathcal{C}_{iJkL} \neq \mathcal{C}_{kLiJ}
//! ```

#[cfg(test)]
pub mod test;

mod almansi_hamel;

pub use almansi_hamel::AlmansiHamel;

use super::*;

/// Required methods for elastic constitutive models.
pub trait Elastic<'a>
where
    Self: Solid<'a>,
{
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma} = J^{-1}\mathbf{P}\cdot\mathbf{F}^T
    /// ```
    fn calculate_cauchy_stress(
        &self,
        deformation_gradient: &DeformationGradient,
    ) -> Result<CauchyStress, ConstitutiveError> {
        Ok(deformation_gradient
            * self.calculate_second_piola_kirchoff_stress(deformation_gradient)?
            * deformation_gradient.transpose()
            / deformation_gradient.determinant())
    }
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{T}_{ijkL} = \frac{\partial\sigma_{ij}}{\partial F_{kL}} = J^{-1} \mathcal{G}_{MNkL} F_{iM} F_{jN} - \sigma_{ij} F_{kL}^{-T} + \left(\delta_{jk}\sigma_{is} + \delta_{ik}\sigma_{js}\right)F_{sL}^{-T}
    /// ```
    fn calculate_cauchy_tangent_stiffness(
        &self,
        deformation_gradient: &DeformationGradient,
    ) -> Result<CauchyTangentStiffness, ConstitutiveError> {
        let deformation_gradient_inverse_transpose = deformation_gradient.inverse_transpose();
        let cauchy_stress = self.calculate_cauchy_stress(deformation_gradient)?;
        let some_stress = &cauchy_stress * &deformation_gradient_inverse_transpose;
        Ok(self
            .calculate_second_piola_kirchoff_tangent_stiffness(deformation_gradient)?
            .contract_first_second_indices_with_second_indices_of(
                deformation_gradient,
                deformation_gradient,
            )
            / deformation_gradient.determinant()
            - CauchyTangentStiffness::dyad_ij_kl(
                &cauchy_stress,
                &deformation_gradient_inverse_transpose,
            )
            + CauchyTangentStiffness::dyad_il_kj(&some_stress, &IDENTITY)
            + CauchyTangentStiffness::dyad_ik_jl(&IDENTITY, &some_stress))
    }
    /// Calculates and returns the first Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathbf{P} = J\boldsymbol{\sigma}\cdot\mathbf{F}^{-T}
    /// ```
    fn calculate_first_piola_kirchoff_stress(
        &self,
        deformation_gradient: &DeformationGradient,
    ) -> Result<FirstPiolaKirchoffStress, ConstitutiveError> {
        Ok(self.calculate_cauchy_stress(deformation_gradient)?
            * deformation_gradient.inverse_transpose()
            * deformation_gradient.determinant())
    }
    /// Calculates and returns the tangent stiffness associated with the first Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathcal{C}_{iJkL} = \frac{\partial P_{iJ}}{\partial F_{kL}} = J \mathcal{T}_{iskL} F_{sJ}^{-T} + P_{iJ} F_{kL}^{-T} - P_{iL} F_{kJ}^{-T}
    /// ```
    fn calculate_first_piola_kirchoff_tangent_stiffness(
        &self,
        deformation_gradient: &DeformationGradient,
    ) -> Result<FirstPiolaKirchoffTangentStiffness, ConstitutiveError> {
        let deformation_gradient_inverse_transpose = deformation_gradient.inverse_transpose();
        let first_piola_kirchoff_stress =
            self.calculate_first_piola_kirchoff_stress(deformation_gradient)?;
        Ok(self
            .calculate_cauchy_tangent_stiffness(deformation_gradient)?
            .contract_second_index_with_first_index_of(&deformation_gradient_inverse_transpose)
            * deformation_gradient.determinant()
            + FirstPiolaKirchoffTangentStiffness::dyad_ij_kl(
                &first_piola_kirchoff_stress,
                &deformation_gradient_inverse_transpose,
            )
            - FirstPiolaKirchoffTangentStiffness::dyad_il_kj(
                &first_piola_kirchoff_stress,
                &deformation_gradient_inverse_transpose,
            ))
    }
    /// Calculates and returns the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathbf{S} = \mathbf{F}^{-1}\cdot\mathbf{P}
    /// ```
    fn calculate_second_piola_kirchoff_stress(
        &self,
        deformation_gradient: &DeformationGradient,
    ) -> Result<SecondPiolaKirchoffStress, ConstitutiveError> {
        Ok(deformation_gradient.inverse()
            * self.calculate_cauchy_stress(deformation_gradient)?
            * deformation_gradient.inverse_transpose()
            * deformation_gradient.determinant())
    }
    /// Calculates and returns the tangent stiffness associated with the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathcal{G}_{IJkL} = \frac{\partial S_{IJ}}{\partial F_{kL}} = \mathcal{C}_{mJkL}F_{mI}^{-T} - S_{LJ}F_{kI}^{-T} = J \mathcal{T}_{mnkL} F_{mI}^{-T} F_{nJ}^{-T} + S_{IJ} F_{kL}^{-T} - S_{IL} F_{kJ}^{-T} -S_{LJ} F_{kI}^{-T}
    /// ```
    fn calculate_second_piola_kirchoff_tangent_stiffness(
        &self,
        deformation_gradient: &DeformationGradient,
    ) -> Result<SecondPiolaKirchoffTangentStiffness, ConstitutiveError> {
        let deformation_gradient_inverse_transpose = deformation_gradient.inverse_transpose();
        let deformation_gradient_inverse = deformation_gradient_inverse_transpose.transpose();
        let second_piola_kirchoff_stress =
            self.calculate_second_piola_kirchoff_stress(deformation_gradient)?;
        Ok(self
            .calculate_cauchy_tangent_stiffness(deformation_gradient)?
            .contract_first_second_indices_with_second_indices_of(
                &deformation_gradient_inverse,
                &deformation_gradient_inverse,
            )
            * deformation_gradient.determinant()
            + SecondPiolaKirchoffTangentStiffness::dyad_ij_kl(
                &second_piola_kirchoff_stress,
                &deformation_gradient_inverse_transpose,
            )
            - SecondPiolaKirchoffTangentStiffness::dyad_il_kj(
                &second_piola_kirchoff_stress,
                &deformation_gradient_inverse_transpose,
            )
            - SecondPiolaKirchoffTangentStiffness::dyad_ik_jl(
                &deformation_gradient_inverse,
                &second_piola_kirchoff_stress,
            ))
    }
    fn solve_uniaxial_tension(
        &self,
        deformation_gradient_11: &Scalar,
    ) -> Result<Scalar, ConstitutiveError> {
        let mut deformation_gradient = IDENTITY_10;
        deformation_gradient[0][0] = *deformation_gradient_11;
        let mut deformation_gradient_22 = 1.0;
        let mut residual;
        let mut residual_abs = 1.0;
        let mut residual_rel = 1.0;
        let mut steps: usize = 0;
        let steps_maximum: usize = 12345;
        while residual_abs >= ABS_TOL && residual_rel >= REL_TOL {
            if steps > steps_maximum {
                return Err(ConstitutiveError::SolveError);
            } else {
                steps += 1;
            }
            deformation_gradient[1][1] = deformation_gradient_22;
            deformation_gradient[2][2] = deformation_gradient_22;
            residual = self.calculate_cauchy_stress(&deformation_gradient)?[1][1];
            residual_abs = residual.abs();
            residual_rel =
                residual_abs / self.calculate_cauchy_stress(&deformation_gradient)?[0][0].abs();
            deformation_gradient_22 -= residual
                / self.calculate_cauchy_tangent_stiffness(&deformation_gradient)?[1][1][1][1];
        }
        Ok(self.calculate_cauchy_stress(&deformation_gradient)?[0][0])
    }
}
