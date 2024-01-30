#[cfg(test)]
pub mod test;

mod almansi_hamel;

use super::*;

/// Required methods for thermoelastic constitutive models.
pub trait ThermoelasticConstitutiveModel<'a>
where
    Self: ConstitutiveModel<'a>
{
    /// Calculates and returns the Cauchy stress.
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient, temperature: &Temperature) -> CauchyStress
    {
        deformation_gradient*self.calculate_second_piola_kirchoff_stress(deformation_gradient, temperature)*deformation_gradient.transpose()/deformation_gradient.determinant()
    }
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient, temperature: &Temperature) -> CauchyTangentStiffness
    {
        let deformation_gradient_inverse_transpose = deformation_gradient.inverse_transpose();
        let cauchy_stress = self.calculate_cauchy_stress(deformation_gradient, temperature);
        let identity = CauchyStress::identity();
        let some_stress = &cauchy_stress * &deformation_gradient_inverse_transpose;
        self.calculate_second_piola_kirchoff_tangent_stiffness(deformation_gradient, temperature).contract_first_second_indices_with_second_indices_of(deformation_gradient, deformation_gradient)/deformation_gradient.determinant() - CauchyTangentStiffness::dyad_ij_kl(&cauchy_stress, &deformation_gradient_inverse_transpose) + CauchyTangentStiffness::dyad_il_kj(&some_stress, &identity) + CauchyTangentStiffness::dyad_ik_jl(&identity, &some_stress)
    }
    /// Calculates and returns the first Piola-Kirchoff stress.
    fn calculate_first_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient, temperature: &Temperature) -> FirstPiolaKirchoffStress
    {
        self.calculate_cauchy_stress(deformation_gradient, temperature)*deformation_gradient.inverse_transpose()*deformation_gradient.determinant()
    }
    /// Calculates and returns the tangent stiffness associated with the first Piola-Kirchoff stress.
    fn calculate_first_piola_kirchoff_tangent_stiffness(&self, deformation_gradient: &DeformationGradient, temperature: &Temperature) -> FirstPiolaKirchoffTangentStiffness
    {
        let deformation_gradient_inverse_transpose = deformation_gradient.inverse_transpose();
        let first_piola_kirchoff_stress = self.calculate_first_piola_kirchoff_stress(deformation_gradient, temperature);
        self.calculate_cauchy_tangent_stiffness(deformation_gradient, temperature).contract_second_index_with_first_index_of(&deformation_gradient_inverse_transpose)*deformation_gradient.determinant() + FirstPiolaKirchoffTangentStiffness::dyad_ij_kl(&first_piola_kirchoff_stress, &deformation_gradient_inverse_transpose) - FirstPiolaKirchoffTangentStiffness::dyad_il_kj(&first_piola_kirchoff_stress, &deformation_gradient_inverse_transpose)
    }
    /// Calculates and returns the second Piola-Kirchoff stress.
    fn calculate_second_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient, temperature: &Temperature) -> SecondPiolaKirchoffStress
    {
        deformation_gradient.inverse()*self.calculate_cauchy_stress(deformation_gradient, temperature)*deformation_gradient.inverse_transpose()*deformation_gradient.determinant()
    }
    /// Calculates and returns the tangent stiffness associated with the second Piola-Kirchoff stress.
    fn calculate_second_piola_kirchoff_tangent_stiffness(&self, deformation_gradient: &DeformationGradient, temperature: &Temperature) -> SecondPiolaKirchoffTangentStiffness
    {
        let deformation_gradient_inverse_transpose = deformation_gradient.inverse_transpose();
        let deformation_gradient_inverse = deformation_gradient_inverse_transpose.transpose();
        let second_piola_kirchoff_stress = self.calculate_second_piola_kirchoff_stress(deformation_gradient, temperature);
        self.calculate_cauchy_tangent_stiffness(deformation_gradient, temperature).contract_first_second_indices_with_second_indices_of(&deformation_gradient_inverse, &deformation_gradient_inverse)*deformation_gradient.determinant() + SecondPiolaKirchoffTangentStiffness::dyad_ij_kl(&second_piola_kirchoff_stress, &deformation_gradient_inverse_transpose) - SecondPiolaKirchoffTangentStiffness::dyad_il_kj(&second_piola_kirchoff_stress, &deformation_gradient_inverse_transpose) - SecondPiolaKirchoffTangentStiffness::dyad_ik_jl(&deformation_gradient_inverse, &second_piola_kirchoff_stress)
    }
    /// Returns the bulk modulus.
    fn get_bulk_modulus(&self) -> &Scalar;
    /// Returns the shear modulus.
    fn get_shear_modulus(&self) -> &Scalar;
    /// Returns the coefficient of thermal expansion.
    fn get_coefficient_of_thermal_expansion(&self) -> &Scalar;
    /// Returns the reference temperature.
    fn get_reference_temperature(&self) -> &Scalar;
}
