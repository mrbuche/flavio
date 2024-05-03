#[cfg(test)]
pub mod test;

use super::*;

pub struct SmithFerrante<'a>
{
    parameters: Parameters<'a>
}

impl<'a> SmithFerrante<'a>
{
    fn get_characteristic_displacement(&self) -> &Scalar
    {
        &self.parameters[0]
    }
    fn get_maximum_normal_traction(&self) -> &Scalar
    {
        &self.parameters[1]
    }
    fn get_weight(&self) -> &Scalar
    {
        &self.parameters[2]
    }
}

impl<'a> Constitutive<'a> for SmithFerrante<'a>
{
    fn new(parameters: Parameters<'a>) -> Self
    {
        Self
        {
            parameters
        }
    }
}

impl<'a> Cohesive<'a> for SmithFerrante<'a>
{
    fn calculate_potential(&self, displacement: &Displacement) -> Scalar
    {
        let characteristic_displacement = self.get_characteristic_displacement();
        let nondimensional_displacement_magnitude = displacement.norm() / characteristic_displacement;
        characteristic_displacement * self.get_maximum_normal_traction() * 1.0_f64.exp() * (1.0 - (1.0 + nondimensional_displacement_magnitude) * (-nondimensional_displacement_magnitude).exp())
        // why do you need this if you cant calculate a useful energy in FEA?
        // also you cant minimize anything due to loading/unloading
        // so in summary, the potential seems useless
        //
        // or can you connect to FEA somehow? try some math on paper
    }
    fn calculate_traction(&self, displacement: &Displacement, normal: &Normal) -> Traction
    {
        let normal_displacement = normal * (displacement * normal);
        let surface_displacement = displacement - normal * (displacement * normal);
        let characteristic_displacement = self.get_characteristic_displacement();
        (surface_displacement * self.get_weight().powi(2) + normal_displacement) * (self.get_maximum_normal_traction() / characteristic_displacement) * (1.0 - displacement.norm() / characteristic_displacement).exp()
    }
}