#[cfg(test)]
mod test;

use super::*;

pub struct SaintVenantKirchoffModel<'a>
{
    parameters: ConstitutiveModelParameters<'a>
}

impl<'a> ConstitutiveModel<'a> for SaintVenantKirchoffModel<'a>
{
    fn new(parameters: ConstitutiveModelParameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

impl<'a> ThermoelasticConstitutiveModel<'a> for SaintVenantKirchoffModel<'a>
{
    fn calculate_second_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient, temperature: &Temperature) -> SecondPiolaKirchoffStress
    {
        let (deviatoric_strain, strain_trace) = ((self.calculate_right_cauchy_green_deformation(deformation_gradient) - RightCauchyGreenDeformation::identity())*0.5).deviatoric_and_trace();
        deviatoric_strain*(2.0*self.get_shear_modulus()) + RightCauchyGreenDeformation::identity()*(self.get_bulk_modulus()*(strain_trace - 3.0*self.get_coefficient_of_thermal_expansion()*(temperature - self.get_reference_temperature())))
    }
    fn calculate_second_piola_kirchoff_tangent_stiffness(&self, deformation_gradient: &DeformationGradient, _: &Temperature) -> SecondPiolaKirchoffTangentStiffness
    {
        let identity = SecondPiolaKirchoffStress::identity();
        let scaled_deformation_gradient_transpose = deformation_gradient.transpose()*self.get_shear_modulus();
        SecondPiolaKirchoffTangentStiffness::dyad_ik_jl(&scaled_deformation_gradient_transpose, &identity) + SecondPiolaKirchoffTangentStiffness::dyad_il_jk(&identity, &scaled_deformation_gradient_transpose) + SecondPiolaKirchoffTangentStiffness::dyad_ij_kl(&(identity*(self.get_bulk_modulus() - 2.0/3.0*self.get_shear_modulus())), deformation_gradient)
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

impl<'a> ThermohyperelasticConstitutiveModel<'a> for SaintVenantKirchoffModel<'a>
{
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient, temperature: &Temperature) -> Scalar
    {
        let strain = (self.calculate_right_cauchy_green_deformation(deformation_gradient) - RightCauchyGreenDeformation::identity())*0.5;
        let strain_trace = strain.trace();
        self.get_shear_modulus()*strain.squared_trace() + 0.5*(self.get_bulk_modulus() - 2.0/3.0*self.get_shear_modulus())*strain_trace.powi(2) - 3.0*self.get_bulk_modulus()*self.get_coefficient_of_thermal_expansion()*(temperature - self.get_reference_temperature())*strain_trace
    }
}
