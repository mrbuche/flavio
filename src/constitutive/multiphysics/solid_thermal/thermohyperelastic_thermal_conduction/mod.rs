//! Thermohyperelastic-thermal conduction constitutive models.

#[cfg(test)]
mod test;

use crate::mechanics::
{
    CauchyStress,
    CauchyTangentStiffness,
    DeformationGradient,
    FirstPiolaKirchoffStress,
    FirstPiolaKirchoffTangentStiffness,
    HeatFlux,
    SecondPiolaKirchoffStress,
    SecondPiolaKirchoffTangentStiffness,
    Scalar,
    TemperatureGradient
};
use super::*;

/// A thermohyperelastic-thermal conduction constitutive model.
#[derive(Debug)]
pub struct ThermohyperelasticThermalConduction<C1, C2>
{
    thermohyperelastic_constitutive_model: C1,
    thermal_conduction_constitutive_model: C2
}

/// Constitutive model implementation of a thermohyperelastic-thermal conduction constitutive model.
impl<'a, C1, C2> Constitutive<'a> for ThermohyperelasticThermalConduction<C1, C2>
{
    /// Dummy method that will panic, use [Self::construct()] instead.
    fn new(_parameters: Parameters<'a>) -> Self
    {
        panic!()
    }
}

/// Solid constitutive model implementation of a thermohyperelastic-thermal conduction constitutive model.
impl<'a, C1, C2> Solid<'a> for ThermohyperelasticThermalConduction<C1, C2>
where
    C1: Thermohyperelastic<'a>,
    C2: ThermalConduction<'a>
{
    fn get_bulk_modulus(&self) -> &Scalar
    {
        self.get_solid_constitutive_model().get_bulk_modulus()
    }
    fn get_shear_modulus(&self) -> &Scalar
    {
        self.get_solid_constitutive_model().get_shear_modulus()
    }
}

/// Thermoelastic constitutive model implementation of a thermohyperelastic-thermal conduction constitutive model.
impl<'a, C1, C2> Thermoelastic<'a> for ThermohyperelasticThermalConduction<C1, C2>
where
    C1: Thermohyperelastic<'a>,
    C2: ThermalConduction<'a>
{
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient, temperature: &Scalar) -> CauchyStress
    {
        self.get_solid_constitutive_model().calculate_cauchy_stress(deformation_gradient, temperature)
    }
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient, temperature: &Scalar) -> CauchyTangentStiffness
    {
        self.get_solid_constitutive_model().calculate_cauchy_tangent_stiffness(deformation_gradient, temperature)
    }
    fn calculate_first_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient, temperature: &Scalar) -> FirstPiolaKirchoffStress
    {
        self.get_solid_constitutive_model().calculate_first_piola_kirchoff_stress(deformation_gradient, temperature)
    }
    fn calculate_first_piola_kirchoff_tangent_stiffness(&self, deformation_gradient: &DeformationGradient, temperature: &Scalar) -> FirstPiolaKirchoffTangentStiffness
    {
        self.get_solid_constitutive_model().calculate_first_piola_kirchoff_tangent_stiffness(deformation_gradient, temperature)
    }
    fn calculate_second_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient, temperature: &Scalar) -> SecondPiolaKirchoffStress
    {
        self.get_solid_constitutive_model().calculate_second_piola_kirchoff_stress(deformation_gradient, temperature)
    }
    fn calculate_second_piola_kirchoff_tangent_stiffness(&self, deformation_gradient: &DeformationGradient, temperature: &Scalar) -> SecondPiolaKirchoffTangentStiffness
    {
        self.get_solid_constitutive_model().calculate_second_piola_kirchoff_tangent_stiffness(deformation_gradient, temperature)
    }
    fn get_coefficient_of_thermal_expansion(&self) -> &Scalar
    {
        self.get_solid_constitutive_model().get_coefficient_of_thermal_expansion()
    }
    fn get_reference_temperature(&self) -> &Scalar
    {
        self.get_solid_constitutive_model().get_reference_temperature()
    }
}

/// Thermohyperelastic constitutive model implementation of a thermohyperelastic-thermal conduction constitutive model.
impl<'a, C1, C2> Thermohyperelastic<'a> for ThermohyperelasticThermalConduction<C1, C2>
where
    C1: Thermohyperelastic<'a>,
    C2: ThermalConduction<'a>
{
    fn calculate_helmholtz_free_energy_density(&'a self, deformation_gradient: &'a DeformationGradient, temperature: &'a Scalar) -> Result<Scalar, ConstitutiveError>
    {
        self.get_solid_constitutive_model().calculate_helmholtz_free_energy_density(deformation_gradient, temperature)
    }
}

/// Thermal constitutive model implementation of a thermohyperelastic-thermal conduction constitutive model.
impl<'a, C1, C2> Thermal<'a> for ThermohyperelasticThermalConduction<C1, C2> {}

/// Thermal conduction constitutive model implementation of a thermhyperoelastic-thermal conduction constitutive model.
impl<'a, C1, C2> ThermalConduction<'a> for ThermohyperelasticThermalConduction<C1, C2>
where
    C1: Thermohyperelastic<'a>,
    C2: ThermalConduction<'a>
{
    fn calculate_heat_flux(&self, temperature_gradient: &TemperatureGradient) -> HeatFlux
    {
        self.get_thermal_constitutive_model().calculate_heat_flux(temperature_gradient)
    }
}

/// Multiphysics constitutive model implementation of a thermohyperelastic-thermal conduction constitutive model.
impl<'a, C1, C2> Multiphysics<'a> for ThermohyperelasticThermalConduction<C1, C2> {}

/// Solid-thermal constitutive model implementation of a thermohyperelastic-thermal conduction constitutive model.
impl<'a, C1, C2> SolidThermal<'a, C1, C2> for ThermohyperelasticThermalConduction<C1, C2>
where
    C1: Thermohyperelastic<'a>,
    C2: ThermalConduction<'a>
{
    fn construct(thermohyperelastic_constitutive_model: C1, thermal_conduction_constitutive_model: C2) -> Self
    {
        Self
        {
            thermohyperelastic_constitutive_model,
            thermal_conduction_constitutive_model
        }
    }
    fn get_solid_constitutive_model(&self) -> &C1
    {
        &self.thermohyperelastic_constitutive_model
    }
    fn get_thermal_constitutive_model(&self) -> &C2
    {
        &self.thermal_conduction_constitutive_model
    }
}
