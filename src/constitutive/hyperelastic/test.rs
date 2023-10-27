macro_rules! test_hyperelastic_constitutive_model
{
    ($hyperelastic_constitutive_model: ident, $constitutive_model_parameters: expr) =>
    {
        use crate::
        {
            constitutive::
            {
                ConstitutiveModel,
                DeformationGradient
            },
            math::TensorRank2Trait
        };
        fn get_hyperelastic_constitutive_model<'a>() -> $hyperelastic_constitutive_model<'a>
        {
            $hyperelastic_constitutive_model::new($constitutive_model_parameters)
        }
        #[test]
        fn get_bulk_modulus()
        {
            assert_eq!(&$constitutive_model_parameters[0], get_hyperelastic_constitutive_model().get_bulk_modulus())
        }
        #[test]
        fn get_shear_modulus()
        {
            assert_eq!( get_hyperelastic_constitutive_model().get_shear_modulus(), &$constitutive_model_parameters[1])
        }
        #[test]
        fn todo()
        {
            todo!("Set up automatic testing, organize properly, include things like objectivity, and test one thing at a time.")
        }
        mod cauchy_stress
        {
            use super::*;
            #[test]
            fn zero()
            {
                get_hyperelastic_constitutive_model().calculate_cauchy_stress(&DeformationGradient::identity()).iter()
                .for_each(|cauchy_stress_i|
                    cauchy_stress_i.iter()
                    .for_each(|cauchy_stress_ij|
                        assert_eq!(cauchy_stress_ij, &0.0)
                    )
                )
            }
        }
        mod helmholtz_free_energy_density
        {
            use super::*;
            #[test]
            fn zero()
            {
                assert_eq!(get_hyperelastic_constitutive_model().calculate_helmholtz_free_energy_density(&DeformationGradient::identity()), 0.0)
            }
        }
    }
}
pub(crate) use test_hyperelastic_constitutive_model;