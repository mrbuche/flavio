//! Thermal constitutive models.

pub mod conduction;

use crate::
{
    mechanics::
    {
        HeatFlux,
        Scalar,
        TemperatureGradient
    }
};

use super::
{
    ConstitutiveModel,
    ConstitutiveModelParameters
};

pub trait ThermalConstitutiveModel<'a>
where
    Self: ConstitutiveModel<'a>
{}
