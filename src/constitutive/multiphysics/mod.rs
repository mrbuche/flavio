//! Multiphysics constitutive models.

use super::
{
    ConstitutiveModel,
    solid::SolidConstitutiveModel,
    thermal::ThermalConstitutiveModel
};

pub trait MultiphysicsConstitutiveModel<'a>
where
    Self: ConstitutiveModel<'a>
{}

pub trait ThermalSolidConstitutiveModel<'a, C1, C2>
where
    C1: SolidConstitutiveModel<'a>,
    C2: ThermalConstitutiveModel<'a>,
    Self: MultiphysicsConstitutiveModel<'a>
{
    fn get_solid_constitutive_model(&self) -> &C1;
    fn get_thermal_constitutive_model(&self) -> &C2;
}
