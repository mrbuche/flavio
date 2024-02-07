//! Viscoelastic constitutive models.

#[cfg(test)]
mod test;

use super::
{
    *, super::fluid::viscous::Viscous
};

/// Required methods for viscoelastic constitutive models.
pub trait Viscoelastic<'a>
where
    Self: Solid<'a> + Viscous<'a>
{
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient, deformation_gradient_dot: &DeformationGradientDot) -> CauchyStress
    {
        deformation_gradient*self.calculate_second_piola_kirchoff_stress(deformation_gradient, deformation_gradient_dot)*deformation_gradient.transpose()/deformation_gradient.determinant()
    }
    fn calculate_first_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient, deformation_gradient_dot: &DeformationGradientDot) -> FirstPiolaKirchoffStress
    {
        self.calculate_cauchy_stress(deformation_gradient, deformation_gradient_dot)*deformation_gradient.inverse_transpose()*deformation_gradient.determinant()
    }
    fn calculate_second_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient, deformation_gradient_dot: &DeformationGradientDot) -> SecondPiolaKirchoffStress
    {
        deformation_gradient.inverse()*self.calculate_cauchy_stress(deformation_gradient, deformation_gradient_dot)*deformation_gradient.inverse_transpose()*deformation_gradient.determinant()
    }
    //
    // appropriate tangent seems to be dP/dFdot, what to call that?
    // or do you need both dP/dF and dP/dFdot?
    //
}