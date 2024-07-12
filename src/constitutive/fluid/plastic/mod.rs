//! Plastic constitutive models.

#[cfg(test)]
pub mod test;

use super::*;

/// Required methods for plastic constitutive models.
pub trait Plastic<'a>
// where
//     Self: Fluid<'a>
{
    fn compute_plastic_stretching_rate(&self, deviatoric_mandel_stress: &MandelStress) -> StretchingRatePlastic
    {
        let deviatoric_mandel_stress_norm = deviatoric_mandel_stress.norm();
        deviatoric_mandel_stress * (self.compute_plastic_flow_rate(&deviatoric_mandel_stress_norm) / deviatoric_mandel_stress_norm)
    }
    fn compute_plastic_flow_rate(&self, deviatoric_mandel_stress_norm: &Scalar) -> Scalar;
    /// Returns the reference plastic flow rate.
    fn get_reference_plastic_flow_rate(&self) -> &Scalar;
    /// Returns the shear strength.
    fn get_shear_strength(&self) -> &Scalar;
}
