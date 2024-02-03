use crate::constitutive::solid::thermohyperelastic::
{
    SaintVenantKirchoffModel,
    test::SAINTVENANTKIRCHOFFPARAMETERS
};
use super::*;

test_thermohyperelastic_thermal_conduction_constitutive_model!(
    SaintVenantKirchoffModel, SAINTVENANTKIRCHOFFPARAMETERS,
    FourierModel, FOURIERPARAMETERS
);

macro_rules! test_thermohyperelastic_thermal_conduction_constitutive_model
{
    ($thermohyperelastic_constitutive_model: ident, $thermohyperelastic_constitutive_model_parameters: expr, 
     $thermal_conduction_constitutive_model: ident, $thermal_conduction_constitutive_model_parameters: expr) =>
    {
        use crate::constitutive::multiphysics::solid_thermal::thermoelastic_thermal_conduction::test::test_thermoelastic_thermal_conduction_constitutive_model;
        test_thermoelastic_thermal_conduction_constitutive_model!(
            ThermohyperelasticThermalConductionConstitutiveModel,
            $thermohyperelastic_constitutive_model, $thermohyperelastic_constitutive_model_parameters, 
            $thermal_conduction_constitutive_model, $thermal_conduction_constitutive_model_parameters
        );
        #[test]
        fn calculate_helmholtz_free_energy_density()
        {
            assert_eq!(
                get_thermoelastic_thermal_conduction_constitutive_model()
                .calculate_helmholtz_free_energy_density(
                    &get_deformation_gradient(), &get_temperature()
                ),
                get_thermoelastic_constitutive_model()
                .calculate_helmholtz_free_energy_density(
                    &get_deformation_gradient(), &get_temperature()
                )
            )
        }
    }
}
pub(crate) use test_thermohyperelastic_thermal_conduction_constitutive_model;