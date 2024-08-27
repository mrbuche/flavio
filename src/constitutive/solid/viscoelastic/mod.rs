//! Viscoelastic constitutive models.
//!
//! ---
//!
//! Viscoelastic constitutive models cannot be defined by a Helmholtz free energy density function and viscous dissipation function.
//! These constitutive models are therefore defined by a relation for the stress as a function of the deformation gradient and rate.
//! Consequently, the rate tangent stiffness associated with the first Piola-Kirchoff stress is not symmetric for viscoelastic models.
//!
//! ```math
//! \mathcal{U}_{iJkL} \neq \mathcal{U}_{kLiJ}
//! ```

#[cfg(test)]
pub mod test;

use super::{super::fluid::viscous::Viscous, *};

/// Required methods for viscoelastic constitutive models.
pub trait Viscoelastic<'a>
where
    Self: Solid<'a> + Viscous<'a>,
{
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma} = J^{-1}\mathbf{P}\cdot\mathbf{F}^T
    /// ```
    fn calculate_cauchy_stress(
        &self,
        deformation_gradient: &DeformationGradient,
        deformation_gradient_rate: &DeformationGradientRate,
    ) -> Result<CauchyStress, ConstitutiveError> {
        Ok(deformation_gradient
            * self.calculate_second_piola_kirchoff_stress(
                deformation_gradient,
                deformation_gradient_rate,
            )?
            * deformation_gradient.transpose()
            / deformation_gradient.determinant())
    }
    /// Calculates and returns the rate tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{V}_{ijkL} = \frac{\partial\sigma_{ij}}{\partial\dot{F}_{kL}} = J^{-1} \mathcal{W}_{MNkL} F_{iM} F_{jN}
    /// ```
    fn calculate_cauchy_rate_tangent_stiffness(
        &self,
        deformation_gradient: &DeformationGradient,
        deformation_gradient_rate: &DeformationGradientRate,
    ) -> Result<CauchyRateTangentStiffness, ConstitutiveError> {
        Ok(self
            .calculate_second_piola_kirchoff_rate_tangent_stiffness(
                deformation_gradient,
                deformation_gradient_rate,
            )?
            .contract_first_second_indices_with_second_indices_of(
                deformation_gradient,
                deformation_gradient,
            )
            / deformation_gradient.determinant())
    }
    /// Calculates and returns the first Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathbf{P} = J\boldsymbol{\sigma}\cdot\mathbf{F}^{-T}
    /// ```
    fn calculate_first_piola_kirchoff_stress(
        &self,
        deformation_gradient: &DeformationGradient,
        deformation_gradient_rate: &DeformationGradientRate,
    ) -> Result<FirstPiolaKirchoffStress, ConstitutiveError> {
        Ok(
            self.calculate_cauchy_stress(deformation_gradient, deformation_gradient_rate)?
                * deformation_gradient.inverse_transpose()
                * deformation_gradient.determinant(),
        )
    }
    /// Calculates and returns the rate tangent stiffness associated with the first Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathcal{U}_{iJkL} = \frac{\partial P_{iJ}}{\partial\dot{F}_{kL}} = J \mathcal{V}_{iskL} F_{sJ}^{-T}
    /// ```
    fn calculate_first_piola_kirchoff_rate_tangent_stiffness(
        &self,
        deformation_gradient: &DeformationGradient,
        deformation_gradient_rate: &DeformationGradientRate,
    ) -> Result<FirstPiolaKirchoffRateTangentStiffness, ConstitutiveError> {
        Ok(self
            .calculate_cauchy_rate_tangent_stiffness(
                deformation_gradient,
                deformation_gradient_rate,
            )?
            .contract_second_index_with_first_index_of(&deformation_gradient.inverse_transpose())
            * deformation_gradient.determinant())
    }
    /// Calculates and returns the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathbf{S} = \mathbf{F}^{-1}\cdot\mathbf{P}
    /// ```
    fn calculate_second_piola_kirchoff_stress(
        &self,
        deformation_gradient: &DeformationGradient,
        deformation_gradient_rate: &DeformationGradientRate,
    ) -> Result<SecondPiolaKirchoffStress, ConstitutiveError> {
        Ok(deformation_gradient.inverse()
            * self.calculate_cauchy_stress(deformation_gradient, deformation_gradient_rate)?
            * deformation_gradient.inverse_transpose()
            * deformation_gradient.determinant())
    }
    /// Calculates and returns the rate tangent stiffness associated with the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathcal{W}_{IJkL} = \frac{\partial S_{IJ}}{\partial\dot{F}_{kL}} = \mathcal{U}_{mJkL}F_{mI}^{-T} = J \mathcal{V}_{mnkL} F_{mI}^{-T} F_{nJ}^{-T}
    /// ```
    fn calculate_second_piola_kirchoff_rate_tangent_stiffness(
        &self,
        deformation_gradient: &DeformationGradient,
        deformation_gradient_rate: &DeformationGradientRate,
    ) -> Result<SecondPiolaKirchoffRateTangentStiffness, ConstitutiveError> {
        let deformation_gradient_inverse = deformation_gradient.inverse();
        Ok(self
            .calculate_cauchy_rate_tangent_stiffness(
                deformation_gradient,
                deformation_gradient_rate,
            )?
            .contract_first_second_indices_with_second_indices_of(
                &deformation_gradient_inverse,
                &deformation_gradient_inverse,
            )
            * deformation_gradient.determinant())
    }
    // COULD GIVE A CONSTANT STRAIN RATE AND A FINAL DEFORMATION GRADIENT AND **IMPLICITLY** INTEGRATE THIS IN TIME FOR THE STRESS HISTORY?
    //
    // a start with some F_n, Fdot_n = f_n, and guess F_{n+1} = F_n
    //   1 solve for a new guess for f_{n+1} = f(F_{n+1}) using sigma_22 = 0 and symmetry
    //   2 solve for a new guess for F_{n+1} using backward Euler and Newton's method
    // b repeatedly alternate 1 and 2 until converged
    // repeat a and b until get to specified F
    //
    // how to chose the timestep?
    // see adaptive RK for an example
    // need to estimate the error over the step and cut back timestep if it is too large
    //
    // should you just use RK4 instead?
    // would also be good to distill into a math function or macro
    fn solve_uniaxial(
        &self,
        deformation_gradient: &DeformationGradient,
        deformation_gradient_rate_11: &Scalar,
    ) -> Result<(DeformationGradientRate, CauchyStress), ConstitutiveError> {
        let mut cauchy_stress = ZERO;
        let mut deformation_gradient_rate = IDENTITY_10;
        deformation_gradient_rate[0][0] = *deformation_gradient_rate_11;
        deformation_gradient_rate[1][1] =
            -deformation_gradient_rate_11 / deformation_gradient[0][0].powf(1.5);
        deformation_gradient_rate[2][2] = deformation_gradient_rate[1][1];
        let mut residual = 0.0;
        let mut residual_abs = 1.0;
        let mut residual_rel = 1.0;
        let mut steps: usize = 0;
        while residual_abs >= ABS_TOL && residual_rel >= REL_TOL {
            if steps > MAXIMUM_STEPS {
                return Err(ConstitutiveError::SolveError);
            } else {
                deformation_gradient_rate[1][1] -= residual
                    / self.calculate_cauchy_rate_tangent_stiffness(
                        deformation_gradient,
                        &deformation_gradient_rate,
                    )?[1][1][1][1];
                deformation_gradient_rate[2][2] = deformation_gradient_rate[1][1];
                cauchy_stress =
                    self.calculate_cauchy_stress(deformation_gradient, &deformation_gradient_rate)?;
                residual = cauchy_stress[1][1];
                residual_abs = residual.abs();
                residual_rel = residual_abs / cauchy_stress[0][0].abs();
                steps += 1;
            }
        }
        Ok((deformation_gradient_rate, cauchy_stress))
    }
}
