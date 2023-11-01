#[cfg(test)]
mod test;

use super::*;

pub struct GentModel<'a>
{
    parameters: ConstitutiveModelParameters<'a>
}

impl<'a> GentModel<'a>
{
    fn get_extensibility(&self) -> &Scalar
    {
        &self.parameters[2]
    }
}

impl<'a> ConstitutiveModel<'a> for GentModel<'a>
{
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        let jacobian = deformation_gradient.determinant();
        let isochoric_left_cauchy_green_deformation = self.calculate_left_cauchy_green_deformation(deformation_gradient)/jacobian.powf(2.0/3.0);
        let (deviatoric_isochoric_left_cauchy_green_deformation, isochoric_left_cauchy_green_deformation_trace) = isochoric_left_cauchy_green_deformation.deviatoric_and_trace();
        let denominator = self.get_extensibility() - isochoric_left_cauchy_green_deformation_trace + 3.0;
        if denominator < 0.0
        {
            panic!()
        }
        (deviatoric_isochoric_left_cauchy_green_deformation*self.get_shear_modulus()*self.get_extensibility()/jacobian)/denominator + LeftCauchyGreenDeformation::identity()*self.get_bulk_modulus()*0.5*(jacobian - 1.0/jacobian)
    }
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        let identity = CauchyStress::identity();
        let (inverse_transpose_deformation_gradient, jacobian) = deformation_gradient.inverse_transpose_and_determinant();
        let isochoric_left_cauchy_green_deformation = self.calculate_left_cauchy_green_deformation(deformation_gradient)/jacobian.powf(2.0/3.0);
        let (deviatoric_isochoric_left_cauchy_green_deformation, isochoric_left_cauchy_green_deformation_trace) = isochoric_left_cauchy_green_deformation.deviatoric_and_trace();
        let denominator = self.get_extensibility() - isochoric_left_cauchy_green_deformation_trace + 3.0;
        if denominator < 0.0
        {
            panic!()
        }
        let prefactor = self.get_shear_modulus()*self.get_extensibility()/jacobian/denominator;
        (CauchyTangentStiffness::dyad_ik_jl(&identity, deformation_gradient) + CauchyTangentStiffness::dyad_il_jk(deformation_gradient, &identity) - CauchyTangentStiffness::dyad_ij_kl(&identity, deformation_gradient)*(2.0/3.0) + CauchyTangentStiffness::dyad_ij_kl(&deviatoric_isochoric_left_cauchy_green_deformation, deformation_gradient)*(2.0/denominator))*(prefactor/jacobian.powf(2.0/3.0)) + CauchyTangentStiffness::dyad_ij_kl(&(identity*(0.5*self.get_bulk_modulus()*(jacobian + 1.0/jacobian)) - deviatoric_isochoric_left_cauchy_green_deformation*prefactor*((5.0 + 2.0*isochoric_left_cauchy_green_deformation_trace/denominator)/3.0)), &inverse_transpose_deformation_gradient)
    }
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar
    {
        let jacobian = deformation_gradient.determinant();
        let factor = (self.calculate_left_cauchy_green_deformation(deformation_gradient).trace()/jacobian.powf(2.0/3.0) - 3.0)/self.get_extensibility();
        if factor > 1.0
        {
            panic!()
        }
        0.5*(-self.get_shear_modulus()*self.get_extensibility()*(1.0 - factor).ln() + self.get_bulk_modulus()*(0.5*(jacobian.powi(2) - 1.0) - jacobian.ln()))
    }
    fn new(parameters: ConstitutiveModelParameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

impl<'a> HyperelasticConstitutiveModel for GentModel<'a>
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
