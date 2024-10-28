//! Hyperelastic constitutive models.
//!
//! ---
//!
//! Hyperelastic constitutive models are completely defined by a Helmholtz free energy density function of the deformation gradient.
//!
//! ```math
//! \mathbf{P}:\dot{\mathbf{F}} - \dot{a}(\mathbf{F}) \geq 0
//! ```
//! Satisfying the second law of thermodynamics (here, equivalent to extremized or zero dissipation) yields a relation for the stress.
//!
//! ```math
//! \mathbf{P} = \frac{\partial a}{\partial\mathbf{F}}
//! ```
//! Consequently, the tangent stiffness associated with the first Piola-Kirchoff stress is symmetric for hyperelastic constitutive models.
//!
//! ```math
//! \mathcal{C}_{iJkL} = \mathcal{C}_{kLiJ}
//! ```

#[cfg(test)]
pub mod test;

mod arruda_boyce;
mod fung;
mod gent;
mod mooney_rivlin;
mod neo_hookean;
mod saint_venant_kirchoff;
mod yeoh;

pub use self::{
    arruda_boyce::ArrudaBoyce, fung::Fung, gent::Gent, mooney_rivlin::MooneyRivlin,
    neo_hookean::NeoHookean, saint_venant_kirchoff::SaintVenantKirchoff, yeoh::Yeoh,
};
use super::{elastic::Elastic, *};
use std::fmt::Debug;

/// Required methods for hyperelastic constitutive models.
pub trait Hyperelastic<'a>
where
    Self: Elastic<'a> + Debug,
{
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a = a(\mathbf{F})
    /// ```
    fn calculate_helmholtz_free_energy_density(
        &self,
        deformation_gradient: &DeformationGradient,
    ) -> Result<Scalar, ConstitutiveError>;
    /// ???
    fn solve_biaxial(
        &self,
        deformation_gradient_11: &Scalar,
        deformation_gradient_22: &Scalar,
    ) -> Result<(DeformationGradient, CauchyStress), ConstitutiveError> {
        let mut cauchy_stress;
        let mut deformation_gradient = IDENTITY_10;
        deformation_gradient[0][0] = *deformation_gradient_11;
        deformation_gradient[1][1] = *deformation_gradient_22;
        deformation_gradient[2][2] = 1.0 / deformation_gradient_11 / deformation_gradient_22;
        let mut residual;
        let mut residual_abs;
        let mut tangent;
        for _ in 0..MAXIMUM_STEPS {
            cauchy_stress = self.calculate_cauchy_stress(&deformation_gradient)?;
            residual = cauchy_stress[2][2];
            residual_abs = residual.abs();
            tangent = self.calculate_cauchy_tangent_stiffness(&deformation_gradient)?[2][2][2][2];
            if residual_abs < ABS_TOL
                || residual_abs / (cauchy_stress[0][0].powi(2) + cauchy_stress[1][1].powi(2)).sqrt()
                    < REL_TOL
            {
                if tangent > 0.0 {
                    return Ok((deformation_gradient, cauchy_stress));
                } else {
                    return Err(ConstitutiveError::NotMinimum(
                        deformation_gradient.copy(),
                        format!("{:?}", &self),
                    ));
                }
            } else {
                deformation_gradient[2][2] -= residual / tangent;
            }
        }
        Err(ConstitutiveError::MaximumStepsReached(
            MAXIMUM_STEPS,
            format!("{:?}", &self),
        ))
    }
    /// ???
    fn solve_uniaxial(
        &self,
        deformation_gradient_11: &Scalar,
    ) -> Result<(DeformationGradient, CauchyStress), ConstitutiveError> {
        let mut cauchy_stress;
        let mut deformation_gradient = IDENTITY_10;
        deformation_gradient[0][0] = *deformation_gradient_11;
        deformation_gradient[1][1] = 1.0 / deformation_gradient_11.sqrt();
        deformation_gradient[2][2] = deformation_gradient[1][1];
        let mut residual;
        let mut residual_abs;
        let mut tangent;
        for _ in 0..MAXIMUM_STEPS {
            cauchy_stress = self.calculate_cauchy_stress(&deformation_gradient)?;
            residual = cauchy_stress[1][1];
            residual_abs = residual.abs();
            tangent = self.calculate_cauchy_tangent_stiffness(&deformation_gradient)?[1][1][1][1];
            if residual_abs < ABS_TOL || residual_abs / cauchy_stress[0][0].abs() < REL_TOL {
                if tangent > 0.0 {
                    return Ok((deformation_gradient, cauchy_stress));
                } else {
                    return Err(ConstitutiveError::NotMinimum(
                        deformation_gradient.copy(),
                        format!("{:?}", &self),
                    ));
                }
            } else {
                deformation_gradient[1][1] -= residual / tangent;
                deformation_gradient[2][2] = deformation_gradient[1][1];
            }
        }
        Err(ConstitutiveError::MaximumStepsReached(
            MAXIMUM_STEPS,
            format!("{:?}", &self),
        ))
    }
}
