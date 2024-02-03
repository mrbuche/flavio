//! Multiphysics constitutive models.

mod solid_thermal;

pub use solid_thermal::
{
    SolidThermalConstitutiveModel,
    thermoelastic_thermal_conduction::ThermoelasticThermalConductionConstitutiveModel,
    thermohyperelastic_thermal_conduction::ThermohyperelasticThermalConductionConstitutiveModel
};

use super::
{
    ConstitutiveModel,
    ConstitutiveModelParameters,
    solid::
    {
        SolidConstitutiveModel,
        thermoelastic::ThermoelasticConstitutiveModel,
        thermohyperelastic::ThermohyperelasticConstitutiveModel
    },
    thermal::
    {
        ThermalConstitutiveModel,
        conduction::ThermalConductionConstitutiveModel
    }
};

/// Required methods for multiphysics constitutive models.
pub trait MultiphysicsConstitutiveModel<'a>
where
    Self: ConstitutiveModel<'a>
{}