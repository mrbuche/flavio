//! Elastic-hyperviscous constitutive models.
//!
//! ---
//!
//! Elastic-hyperviscous constitutive models are defined by an elastic stress tensor function and a viscous dissipation function.
//!
//! ```math
//! \mathbf{P}:\dot{\mathbf{F}} - \mathbf{P}^e(\mathbf{F}):\dot{\mathbf{F}} - \phi(\mathbf{F},\dot{\mathbf{F}}) \geq 0
//! ```
//! Satisfying the second law of thermodynamics though a minimum viscous dissipation principal yields a relation for the stress.
//!
//! ```math
//! \mathbf{P} = \mathbf{P}^e + \frac{\partial\phi}{\partial\dot{\mathbf{F}}}
//! ```
//! Consequently, the rate tangent stiffness associated with the first Piola-Kirchoff stress is symmetric for elastic-hyperviscous models.
//!
//! ```math
//! \mathcal{U}_{iJkL} = \mathcal{U}_{kLiJ}
//! ```

#[cfg(test)]
pub mod test;

mod almansi_hamel;

pub use almansi_hamel::AlmansiHamel;

use super::{super::fluid::viscous::Viscous, viscoelastic::Viscoelastic, *};
use std::fmt::Debug;

/// Required methods for elastic-hyperviscous constitutive models.
pub trait ElasticHyperviscous<'a>
where
    Self: Viscoelastic<'a> + Debug,
{
    /// Calculates and returns the dissipation potential.
    ///
    /// ```math
    /// \mathbf{P}^e(\mathbf{F}):\dot{\mathbf{F}} + \phi(\mathbf{F},\dot{\mathbf{F}})
    /// ```
    fn calculate_dissipation_potential(
        &self,
        deformation_gradient: &DeformationGradient,
        deformation_gradient_rate: &DeformationGradientRate,
    ) -> Result<Scalar, ConstitutiveError> {
        Ok(self
            .calculate_first_piola_kirchoff_stress(deformation_gradient, &ZERO_10)?
            .full_contraction(deformation_gradient_rate)
            + self
                .calculate_viscous_dissipation(deformation_gradient, deformation_gradient_rate)?)
    }
    /// Calculates and returns the viscous dissipation.
    ///
    /// ```math
    /// \phi = \phi(\mathbf{F},\dot{\mathbf{F}})
    /// ```
    fn calculate_viscous_dissipation(
        &self,
        deformation_gradient: &DeformationGradient,
        deformation_gradient_rate: &DeformationGradientRate,
    ) -> Result<Scalar, ConstitutiveError>;
    /// Solve for the unknown components of the Cauchy stress and deformation gradient under uniaxial stress.
    fn solve_uniaxial<const W: usize>(
        &self,
        deformation_gradient_rate_11: impl Fn(Scalar) -> Scalar,
        evaluation_times: [Scalar; W],
    ) -> Result<(DeformationGradients<W>, CauchyStresses<W>), ConstitutiveError> {
        let mut cauchy_stresses = CauchyStresses::<W>::zero();
        let mut deformation_gradients = DeformationGradients::<W>::identity();
        let time_steps = evaluation_times.windows(2).map(|time| time[1] - time[0]);
        for ((index, time_step), time) in time_steps.enumerate().zip(evaluation_times.into_iter()) {
            (deformation_gradients[index + 1], cauchy_stresses[index + 1]) = self
                .solve_uniaxial_inner(
                    &deformation_gradients[index],
                    deformation_gradient_rate_11(time),
                    time_step,
                )?;
        }
        Ok((deformation_gradients, cauchy_stresses))
    }
    #[doc(hidden)]
    fn solve_uniaxial_inner(
        &self,
        deformation_gradient_previous: &DeformationGradient,
        deformation_gradient_rate_11: Scalar,
        time_step: Scalar,
    ) -> Result<(DeformationGradient, CauchyStress), ConstitutiveError> {
        let mut cauchy_stress;
        let mut deformation_gradient = IDENTITY_10;
        let mut deformation_gradient_rate;
        let mut residual;
        let mut tangent;
        for _ in 0..MAXIMUM_STEPS {
            (deformation_gradient_rate, cauchy_stress) = self
                .solve_uniaxial_inner_inner(&deformation_gradient, &deformation_gradient_rate_11)?;
            residual = deformation_gradient.copy()
                - deformation_gradient_previous
                - &deformation_gradient_rate * time_step;
            tangent = IDENTITY_1010
                - TensorRank4::dyad_ik_jl(
                    &(&deformation_gradient_rate * deformation_gradient.inverse()),
                    &IDENTITY_00,
                ) * time_step;
            if residual.norm() < ABS_TOL {
                return Ok((deformation_gradient, cauchy_stress));
            } else {
                deformation_gradient -= residual / tangent;
            }
        }
        Err(ConstitutiveError::MaximumStepsReached(
            MAXIMUM_STEPS,
            format!("{:?}", &self),
        ))
    }
    #[doc(hidden)]
    fn solve_uniaxial_inner_inner(
        &self,
        deformation_gradient: &DeformationGradient,
        deformation_gradient_rate_11: &Scalar,
    ) -> Result<(DeformationGradientRate, CauchyStress), ConstitutiveError> {
        let mut cauchy_stress;
        let mut deformation_gradient_rate = ZERO_10;
        deformation_gradient_rate[0][0] = *deformation_gradient_rate_11;
        deformation_gradient_rate[1][1] =
            -deformation_gradient_rate_11 / deformation_gradient[0][0].powf(1.5);
        deformation_gradient_rate[2][2] = deformation_gradient_rate[1][1];
        let mut residual;
        let mut residual_abs;
        let mut tangent;
        for _ in 0..MAXIMUM_STEPS {
            cauchy_stress =
                self.calculate_cauchy_stress(deformation_gradient, &deformation_gradient_rate)?;
            residual = cauchy_stress[1][1];
            residual_abs = residual.abs();
            tangent = self.calculate_cauchy_rate_tangent_stiffness(
                deformation_gradient,
                &deformation_gradient_rate,
            )?[1][1][1][1];
            if residual_abs < ABS_TOL || residual_abs / cauchy_stress[0][0].abs() < REL_TOL {
                if tangent > 0.0 {
                    return Ok((deformation_gradient_rate, cauchy_stress));
                } else {
                    return Err(ConstitutiveError::NotMinimum(
                        format!(
                            "From deformation gradient rate: {}.\n\
                                 For deformation gradient: {}.",
                            deformation_gradient_rate, deformation_gradient
                        ),
                        format!("{:?}", &self),
                    ));
                }
            } else {
                deformation_gradient_rate[1][1] -= residual / tangent;
                deformation_gradient_rate[2][2] = deformation_gradient_rate[1][1];
            }
        }
        Err(ConstitutiveError::MaximumStepsReached(
            MAXIMUM_STEPS,
            format!("{:?}", &self),
        ))
    }
}
