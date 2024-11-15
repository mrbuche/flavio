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
use crate::math::optimize::{NewtonRaphson, SecondOrder};
use std::fmt::Debug;

/// Possible applied loads.
pub enum AppliedLoad {
    /// Uniaxial stress given $`F_{11}`$.
    UniaxialStress(Scalar),
    /// Biaxial stress given $`F_{11}`$ and $`F_{22}`$.
    BiaxialStress(Scalar, Scalar),
}

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
    /// Solve for the unknown components of the Cauchy stress and deformation gradient under an applied load.
    fn solve(
        &self,
        applied_load: AppliedLoad,
    ) -> Result<(DeformationGradient, CauchyStress), ConstitutiveError> {
        let optimization = NewtonRaphson {
            ..Default::default()
        };
        let deformation_gradient = match applied_load {
            AppliedLoad::UniaxialStress(deformation_gradient_11) => {
                let deformation_gradient_33 = optimization.minimize::<1, 0>(
                    |deformation_gradient_33: &Scalar| {
                        Ok(self.calculate_cauchy_stress(&DeformationGradient::new([
                            [deformation_gradient_11, 0.0, 0.0],
                            [0.0, *deformation_gradient_33, 0.0],
                            [0.0, 0.0, *deformation_gradient_33],
                        ]))?[2][2])
                    },
                    |deformation_gradient_33: &Scalar| {
                        Ok(
                            self.calculate_cauchy_tangent_stiffness(&DeformationGradient::new([
                                [deformation_gradient_11, 0.0, 0.0],
                                [0.0, *deformation_gradient_33, 0.0],
                                [0.0, 0.0, *deformation_gradient_33],
                            ]))?[2][2][2][2],
                        )
                    },
                    1.0 / deformation_gradient_11.sqrt(),
                    None,
                    None,
                )?;
                DeformationGradient::new([
                    [deformation_gradient_11, 0.0, 0.0],
                    [0.0, deformation_gradient_33, 0.0],
                    [0.0, 0.0, deformation_gradient_33],
                ])
            }
            AppliedLoad::BiaxialStress(deformation_gradient_11, deformation_gradient_22) => {
                let deformation_gradient_33 = optimization.minimize::<1, 0>(
                    |deformation_gradient_33: &Scalar| {
                        Ok(self.calculate_cauchy_stress(&DeformationGradient::new([
                            [deformation_gradient_11, 0.0, 0.0],
                            [0.0, deformation_gradient_22, 0.0],
                            [0.0, 0.0, *deformation_gradient_33],
                        ]))?[2][2])
                    },
                    |deformation_gradient_33: &Scalar| {
                        Ok(
                            self.calculate_cauchy_tangent_stiffness(&DeformationGradient::new([
                                [deformation_gradient_11, 0.0, 0.0],
                                [0.0, deformation_gradient_22, 0.0],
                                [0.0, 0.0, *deformation_gradient_33],
                            ]))?[2][2][2][2],
                        )
                    },
                    1.0 / deformation_gradient_11 / deformation_gradient_22,
                    None,
                    None,
                )?;
                DeformationGradient::new([
                    [deformation_gradient_11, 0.0, 0.0],
                    [0.0, deformation_gradient_22, 0.0],
                    [0.0, 0.0, deformation_gradient_33],
                ])
            }
        };
        let cauchy_stress = self.calculate_cauchy_stress(&deformation_gradient)?;
        Ok((deformation_gradient, cauchy_stress))
    }
}
