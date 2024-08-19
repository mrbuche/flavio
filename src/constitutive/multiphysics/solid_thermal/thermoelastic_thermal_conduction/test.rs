use super::*;
use crate::constitutive::solid::thermoelastic::{test::ALMANSIHAMELPARAMETERS, AlmansiHamel};

test_thermoelastic_thermal_conduction_constitutive_model!(
    ThermoelasticThermalConduction,
    AlmansiHamel,
    ALMANSIHAMELPARAMETERS,
    Fourier,
    FOURIERPARAMETERS
);

macro_rules! test_thermoelastic_thermal_conduction_constitutive_model
{
    ($thermoelastic_thermal_conduction_constitutive_model: ident,
     $thermoelastic_constitutive_model: ident, $thermoelastic_constitutive_model_parameters: expr,
     $thermal_conduction_constitutive_model: ident, $thermal_conduction_constitutive_model_parameters: expr) =>
    {
        use crate::
        {
            constitutive::
            {
                thermal::conduction::
                {
                    Fourier,
                    test::FOURIERPARAMETERS
                }
            },
            mechanics::test::
            {
                get_deformation_gradient,
                get_temperature,
                get_temperature_gradient
            }
        };
        fn get_thermoelastic_constitutive_model<'a>() -> $thermoelastic_constitutive_model<'a>
        {
            $thermoelastic_constitutive_model::new($thermoelastic_constitutive_model_parameters)
        }
        fn get_thermal_conduction_constitutive_model<'a>() -> $thermal_conduction_constitutive_model<'a>
        {
            $thermal_conduction_constitutive_model::new($thermal_conduction_constitutive_model_parameters)
        }
        fn get_thermoelastic_thermal_conduction_constitutive_model<'a>() -> $thermoelastic_thermal_conduction_constitutive_model<$thermoelastic_constitutive_model<'a>, $thermal_conduction_constitutive_model<'a>>
        {
            $thermoelastic_thermal_conduction_constitutive_model::construct(
                get_thermoelastic_constitutive_model(),
                get_thermal_conduction_constitutive_model()
            )
        }
        #[test]
        #[should_panic]
        fn new()
        {
            $thermoelastic_thermal_conduction_constitutive_model::<$thermoelastic_constitutive_model, $thermal_conduction_constitutive_model>::new(FOURIERPARAMETERS);
        }
        #[test]
        fn get_bulk_modulus()
        {
            assert_eq!(
                get_thermoelastic_thermal_conduction_constitutive_model().get_bulk_modulus(),
                get_thermoelastic_constitutive_model().get_bulk_modulus()
            )
        }
        #[test]
        fn get_shear_modulus()
        {
            assert_eq!(
                get_thermoelastic_thermal_conduction_constitutive_model().get_shear_modulus(),
                get_thermoelastic_constitutive_model().get_shear_modulus()
            )
        }
        #[test]
        fn get_coefficient_of_thermal_expansion()
        {
            assert_eq!(
                get_thermoelastic_thermal_conduction_constitutive_model().get_coefficient_of_thermal_expansion(),
                get_thermoelastic_constitutive_model().get_coefficient_of_thermal_expansion()
            )
        }
        #[test]
        fn get_reference_temperature()
        {
            assert_eq!(
                get_thermoelastic_thermal_conduction_constitutive_model().get_reference_temperature(),
                get_thermoelastic_constitutive_model().get_reference_temperature()
            )
        }
        #[test]
        fn calculate_cauchy_stress()
        {
            get_thermoelastic_thermal_conduction_constitutive_model()
            .calculate_cauchy_stress(
                &get_deformation_gradient(), &get_temperature()
            ).expect("the unexpected").iter().zip(
                get_thermoelastic_constitutive_model()
                .calculate_cauchy_stress(
                    &get_deformation_gradient(), &get_temperature()
                ).expect("the unexpected").iter()
            ).for_each(|(cauchy_stress_i, cauchy_stress_solid_i)|
                cauchy_stress_i.iter()
                .zip(cauchy_stress_solid_i.iter())
                .for_each(|(cauchy_stress_ij, cauchy_stress_solid_ij)|
                    assert_eq!(cauchy_stress_ij, cauchy_stress_solid_ij)
                )
            )
        }
        #[test]
        fn calculate_cauchy_tangent_stiffness()
        {
            get_thermoelastic_thermal_conduction_constitutive_model()
            .calculate_cauchy_tangent_stiffness(
                &get_deformation_gradient(), &get_temperature()
            ).expect("the unexpected").iter().zip(
                get_thermoelastic_constitutive_model()
                .calculate_cauchy_tangent_stiffness(
                    &get_deformation_gradient(), &get_temperature()
                ).expect("the unexpected").iter()
            ).for_each(|(cauchy_tangent_stiffness_i, cauchy_tangent_stiffness_solid_i)|
                cauchy_tangent_stiffness_i.iter()
                .zip(cauchy_tangent_stiffness_solid_i.iter())
                .for_each(|(cauchy_tangent_stiffness_ij, cauchy_tangent_stiffness_solid_ij)|
                    cauchy_tangent_stiffness_ij.iter()
                    .zip(cauchy_tangent_stiffness_solid_ij.iter())
                    .for_each(|(cauchy_tangent_stiffness_ijk, cauchy_tangent_stiffness_solid_ijk)|
                        cauchy_tangent_stiffness_ijk.iter()
                        .zip(cauchy_tangent_stiffness_solid_ijk.iter())
                        .for_each(|(cauchy_tangent_stiffness_ijkl, cauchy_tangent_stiffness_solid_ijkl)|
                            assert_eq!(cauchy_tangent_stiffness_ijkl, cauchy_tangent_stiffness_solid_ijkl)
                        )
                    )
                )
            )
        }
        #[test]
        fn calculate_first_piola_kirchoff_stress()
        {
            get_thermoelastic_thermal_conduction_constitutive_model()
            .calculate_first_piola_kirchoff_stress(
                &get_deformation_gradient(), &get_temperature()
            ).expect("the unexpected").iter().zip(
                get_thermoelastic_constitutive_model()
                .calculate_first_piola_kirchoff_stress(
                    &get_deformation_gradient(), &get_temperature()
                ).expect("the unexpected").iter()
            ).for_each(|(first_piola_kirchoff_stress_i, first_piola_kirchoff_stress_solid_i)|
                first_piola_kirchoff_stress_i.iter()
                .zip(first_piola_kirchoff_stress_solid_i.iter())
                .for_each(|(first_piola_kirchoff_stress_ij, first_piola_kirchoff_stress_solid_ij)|
                    assert_eq!(first_piola_kirchoff_stress_ij, first_piola_kirchoff_stress_solid_ij)
                )
            )
        }
        #[test]
        fn calculate_first_piola_kirchoff_tangent_stiffness()
        {
            get_thermoelastic_thermal_conduction_constitutive_model()
            .calculate_first_piola_kirchoff_tangent_stiffness(
                &get_deformation_gradient(), &get_temperature()
            ).expect("the unexpected").iter().zip(
                get_thermoelastic_constitutive_model()
                .calculate_first_piola_kirchoff_tangent_stiffness(
                    &get_deformation_gradient(), &get_temperature()
                ).expect("the unexpected").iter()
            ).for_each(|(first_piola_kirchoff_tangent_stiffness_i, first_piola_kirchoff_tangent_stiffness_solid_i)|
                first_piola_kirchoff_tangent_stiffness_i.iter()
                .zip(first_piola_kirchoff_tangent_stiffness_solid_i.iter())
                .for_each(|(first_piola_kirchoff_tangent_stiffness_ij, first_piola_kirchoff_tangent_stiffness_solid_ij)|
                    first_piola_kirchoff_tangent_stiffness_ij.iter()
                    .zip(first_piola_kirchoff_tangent_stiffness_solid_ij.iter())
                    .for_each(|(first_piola_kirchoff_tangent_stiffness_ijk, first_piola_kirchoff_tangent_stiffness_solid_ijk)|
                        first_piola_kirchoff_tangent_stiffness_ijk.iter()
                        .zip(first_piola_kirchoff_tangent_stiffness_solid_ijk.iter())
                        .for_each(|(first_piola_kirchoff_tangent_stiffness_ijkl, first_piola_kirchoff_tangent_stiffness_solid_ijkl)|
                            assert_eq!(first_piola_kirchoff_tangent_stiffness_ijkl, first_piola_kirchoff_tangent_stiffness_solid_ijkl)
                        )
                    )
                )
            )
        }
        #[test]
        fn calculate_heat_flux()
        {
            get_thermoelastic_thermal_conduction_constitutive_model()
            .calculate_heat_flux(
                &get_temperature_gradient()
            ).iter().zip(
                get_thermal_conduction_constitutive_model()
                .calculate_heat_flux(
                    &get_temperature_gradient()
                ).iter()
            ).for_each(|(heat_flux_i, heat_flux_thermal_i)|
                assert_eq!(heat_flux_i, heat_flux_thermal_i)
            )
        }
        #[test]
        fn calculate_second_piola_kirchoff_stress()
        {
            get_thermoelastic_thermal_conduction_constitutive_model()
            .calculate_second_piola_kirchoff_stress(
                &get_deformation_gradient(), &get_temperature()
            ).expect("the unexpected").iter().zip(
                get_thermoelastic_constitutive_model()
                .calculate_second_piola_kirchoff_stress(
                    &get_deformation_gradient(), &get_temperature()
                ).expect("the unexpected").iter()
            ).for_each(|(second_piola_kirchoff_stress_i, second_piola_kirchoff_stress_solid_i)|
                second_piola_kirchoff_stress_i.iter()
                .zip(second_piola_kirchoff_stress_solid_i.iter())
                .for_each(|(second_piola_kirchoff_stress_ij, second_piola_kirchoff_stress_solid_ij)|
                    assert_eq!(second_piola_kirchoff_stress_ij, second_piola_kirchoff_stress_solid_ij)
                )
            )
        }
        #[test]
        fn calculate_second_piola_kirchoff_tangent_stiffness()
        {
            get_thermoelastic_thermal_conduction_constitutive_model()
            .calculate_second_piola_kirchoff_tangent_stiffness(
                &get_deformation_gradient(), &get_temperature()
            ).expect("the unexpected").iter().zip(
                get_thermoelastic_constitutive_model()
                .calculate_second_piola_kirchoff_tangent_stiffness(
                    &get_deformation_gradient(), &get_temperature()
                ).expect("the unexpected").iter()
            ).for_each(|(second_piola_kirchoff_tangent_stiffness_i, second_piola_kirchoff_tangent_stiffness_solid_i)|
                second_piola_kirchoff_tangent_stiffness_i.iter()
                .zip(second_piola_kirchoff_tangent_stiffness_solid_i.iter())
                .for_each(|(second_piola_kirchoff_tangent_stiffness_ij, second_piola_kirchoff_tangent_stiffness_solid_ij)|
                    second_piola_kirchoff_tangent_stiffness_ij.iter()
                    .zip(second_piola_kirchoff_tangent_stiffness_solid_ij.iter())
                    .for_each(|(second_piola_kirchoff_tangent_stiffness_ijk, second_piola_kirchoff_tangent_stiffness_solid_ijk)|
                        second_piola_kirchoff_tangent_stiffness_ijk.iter()
                        .zip(second_piola_kirchoff_tangent_stiffness_solid_ijk.iter())
                        .for_each(|(second_piola_kirchoff_tangent_stiffness_ijkl, second_piola_kirchoff_tangent_stiffness_solid_ijkl)|
                            assert_eq!(second_piola_kirchoff_tangent_stiffness_ijkl, second_piola_kirchoff_tangent_stiffness_solid_ijkl)
                        )
                    )
                )
            )
        }
        #[test]
        fn size()
        {
            assert_eq!(
                std::mem::size_of::<ThermoelasticThermalConduction<
                    $thermoelastic_constitutive_model, $thermal_conduction_constitutive_model
                >>(), 2 * std::mem::size_of::<crate::constitutive::Parameters>()
            )
        }
    }
}
pub(crate) use test_thermoelastic_thermal_conduction_constitutive_model;
