#[cfg(test)]
mod test;

use super::*;

pub struct AlmansiHamelModel<'a>
{
    parameters: ConstitutiveModelParameters<'a>
}

impl<'a> ConstitutiveModel<'a> for AlmansiHamelModel<'a>
{
    fn new(parameters: ConstitutiveModelParameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

impl<'a> ThermoelasticConstitutiveModel for AlmansiHamelModel<'a>
{
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient, temperature: &Temperature) -> CauchyStress
    {
        let identity = LeftCauchyGreenDeformation::identity();
        let (inverse_deformation_gradient, jacobian) = deformation_gradient.inverse_and_determinant();
        let strain = &identity * 1.0 - inverse_deformation_gradient.transpose() * &inverse_deformation_gradient;
        let (deviatoric_strain, strain_trace) = strain.deviatoric_and_trace();
        deviatoric_strain * (self.get_shear_modulus() / jacobian) + identity * (self.get_bulk_modulus() * (strain_trace / 2.0 / jacobian - 3.0*self.get_coefficient_of_thermal_expansion()*(temperature - self.get_reference_temperature())))
    }
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient, _: &Temperature) -> CauchyTangentStiffness
    {
        let identity = LeftCauchyGreenDeformation::identity();
        let (inverse_transpose_deformation_gradient, jacobian) = deformation_gradient.inverse_transpose_and_determinant();
        let inverse_left_cauchy_green_deformation = &inverse_transpose_deformation_gradient * inverse_transpose_deformation_gradient.transpose();
        let strain = &identity * 1.0 - &inverse_left_cauchy_green_deformation;
        let (deviatoric_strain, strain_trace) = strain.deviatoric_and_trace();
        (CauchyTangentStiffness::dyad_il_jk(&inverse_transpose_deformation_gradient, &inverse_left_cauchy_green_deformation) + CauchyTangentStiffness::dyad_ik_jl(&inverse_left_cauchy_green_deformation, &inverse_transpose_deformation_gradient)) * (self.get_shear_modulus() / jacobian)+ CauchyTangentStiffness::dyad_ij_kl(&identity,&(inverse_left_cauchy_green_deformation * &inverse_transpose_deformation_gradient*((self.get_bulk_modulus() - self.get_shear_modulus() * 2.0 / 3.0) / jacobian))) - CauchyTangentStiffness::dyad_ij_kl(&(deviatoric_strain * (self.get_shear_modulus() / jacobian) + identity * (self.get_bulk_modulus() * strain_trace / 2.0 / jacobian)), &inverse_transpose_deformation_gradient)
    }
    fn get_bulk_modulus(&self) -> &Scalar
    {
        &self.parameters[0]
    }
    fn get_shear_modulus(&self) -> &Scalar
    {
        &self.parameters[1]
    }
    fn get_coefficient_of_thermal_expansion(&self) -> &Scalar
    {
        &self.parameters[2]
    }
    fn get_reference_temperature(&self) -> &Scalar
    {
        &self.parameters[3]
    }
}
