//! Multiphysics constitutive models.

mod solid_thermal;

pub use solid_thermal::
{
    SolidThermal,
    thermoelastic_thermal_conduction::ThermoelasticThermalConduction,
    thermohyperelastic_thermal_conduction::ThermohyperelasticThermalConduction
};

use super::
{
    Constitutive,
    ConstitutiveError,
    Parameters,
    solid::
    {
        Solid,
        thermoelastic::Thermoelastic,
        thermohyperelastic::Thermohyperelastic
    },
    thermal::
    {
        Thermal,
        conduction::ThermalConduction
    }
};

/// Required methods for multiphysics constitutive models.
pub trait Multiphysics<'a>
where
    Self: Constitutive<'a>
{}
