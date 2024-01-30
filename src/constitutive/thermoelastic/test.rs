use crate::
{
    constitutive::test::NEOHOOKEANPARAMETERS,
    mechanics::Scalar
};
pub const ALMANSIHAMELPARAMETERS: &[Scalar; 4] = &[NEOHOOKEANPARAMETERS[0], NEOHOOKEANPARAMETERS[1], 1.0, 1.0];

macro_rules! test_thermoelastic_constitutive_model
{
    ($thermoelastic_constitutive_model: ident, $thermoelastic_constitutive_model_parameters: expr, $thermoelastic_constitutive_model_constructed: expr) =>
    {
        fn get_thermoelastic_constitutive_model<'a>() -> $thermoelastic_constitutive_model<'a>
        {
            $thermoelastic_constitutive_model::new($thermoelastic_constitutive_model_parameters)
        }
        #[test]
        fn get_coefficient_of_thermal_expansion()
        {
            assert_eq!(
                &$thermoelastic_constitutive_model_parameters[2],
                get_thermoelastic_constitutive_model().get_coefficient_of_thermal_expansion()
            )
        }
        #[test]
        fn get_reference_temperature()
        {
            assert_eq!(
                &$thermoelastic_constitutive_model_parameters[3],
                get_thermoelastic_constitutive_model().get_reference_temperature()
            )
        }
        #[test]
        fn coefficient_of_thermal_expansion()
        {
            let model = get_thermoelastic_constitutive_model();
            let deformation_gradient = DeformationGradient::identity();
            let temperature = model.get_reference_temperature() - EPSILON;
            let first_piola_kirchoff_stress = model.calculate_first_piola_kirchoff_stress(&deformation_gradient, &temperature);
            let compare = 3.0*model.get_bulk_modulus()*EPSILON;
            assert!((first_piola_kirchoff_stress[0][0]/compare - model.get_coefficient_of_thermal_expansion()).abs() < EPSILON);
            assert!((first_piola_kirchoff_stress[1][1]/compare - model.get_coefficient_of_thermal_expansion()).abs() < EPSILON);
            assert!((first_piola_kirchoff_stress[2][2]/compare - model.get_coefficient_of_thermal_expansion()).abs() < EPSILON);
            assert_eq!(first_piola_kirchoff_stress[0][1], 0.0);
            assert_eq!(first_piola_kirchoff_stress[0][2], 0.0);
            assert_eq!(first_piola_kirchoff_stress[1][0], 0.0);
            assert_eq!(first_piola_kirchoff_stress[1][2], 0.0);
            assert_eq!(first_piola_kirchoff_stress[2][0], 0.0);
            assert_eq!(first_piola_kirchoff_stress[2][1], 0.0);
        }
        #[test]
        fn reference_temperature()
        {
            todo!("Maybe use this in deformed/undeformed to check same as elastic counterpart at reference temperature?")
        }
        #[test]
        fn size()
        {
            assert_eq!(
                std::mem::size_of::<$thermoelastic_constitutive_model>(),
                std::mem::size_of::<crate::constitutive::ConstitutiveModelParameters>()
            )
        }
        crate::constitutive::thermoelastic::test::test_thermoelastic_constitutive_model_constructed!($thermoelastic_constitutive_model_constructed);
    }
}
pub(crate) use test_thermoelastic_constitutive_model;
macro_rules! test_thermoelastic_constitutive_model_constructed
{
    ($thermoelastic_constitutive_model_constructed: expr) =>
    {
        use crate::
        {
            EPSILON,
            constitutive::
            {
                ConstitutiveModel
            },
            math::
            {
                TensorRank2Trait
            },
            mechanics::
            {
                DeformationGradient
            }
        };
        #[test]
        fn figure_out_how_to_reuse_elastic_hyperelastic_tests()
        {
            todo!()
        }
    }
}
pub(crate) use test_thermoelastic_constitutive_model_constructed;
macro_rules! test_thermoelastic_only_constitutive_model_constructed
{
    ($thermoelastic_constitutive_model_constructed: expr) =>
    {
        #[test]
        fn todo()
        {
            todo!()
        }
    }
}
pub(crate) use test_thermoelastic_only_constitutive_model_constructed;