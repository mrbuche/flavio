macro_rules! test_hyperelastic_constitutive_model
{
    ($hyperelastic_constitutive_model: ident, $constitutive_model_parameters: expr) =>
    {
        const EPSILON: Scalar = 1e-3;
        use crate::
        {
            constitutive::
            {
                ConstitutiveModel,
                DeformationGradient,
                Scalar
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
            mod deformed
            {
                use super::*;

            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn bulk_modulus()
                {
                    let deformation_gradient = DeformationGradient::identity() * (1.0 + EPSILON).powf(1.0/3.0);
                    let model = get_hyperelastic_constitutive_model();
                    let cauchy_stress = model.calculate_cauchy_stress(&deformation_gradient);
                    assert!((3.0*EPSILON*model.get_bulk_modulus()/cauchy_stress.trace() - 1.0).abs() < EPSILON)
                }
                #[test]
                fn shear_modulus()
                {
                    let mut deformation_gradient = DeformationGradient::identity();
                    deformation_gradient[0][1] = EPSILON;
                    let model = get_hyperelastic_constitutive_model();
                    let cauchy_stress = model.calculate_cauchy_stress(&deformation_gradient);
                    assert!((EPSILON*model.get_shear_modulus()/cauchy_stress[0][1] - 1.0).abs() < EPSILON)
                }
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
        }
        mod helmholtz_free_energy_density
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn positive()
                {
                    assert!(get_hyperelastic_constitutive_model().calculate_helmholtz_free_energy_density(&get_deformation_gradient()) > 0.0)
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn zero()
                {
                    assert_eq!(get_hyperelastic_constitutive_model().calculate_helmholtz_free_energy_density(&DeformationGradient::identity()), 0.0)
                }
            }
        }
        mod first_piola_kirchoff_stress
        {
            use super::*;
            mod deformed
            {
                use super::*;

            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn bulk_modulus()
                {
                    let deformation_gradient = DeformationGradient::identity() * (1.0 + EPSILON).powf(1.0/3.0);
                    let model = get_hyperelastic_constitutive_model();
                    let first_piola_kirchoff_stress = model.calculate_first_piola_kirchoff_stress(&deformation_gradient);
                    assert!((3.0*EPSILON*model.get_bulk_modulus()/first_piola_kirchoff_stress.trace() - 1.0).abs() < EPSILON)
                }
                #[test]
                fn shear_modulus()
                {
                    let mut deformation_gradient = DeformationGradient::identity();
                    deformation_gradient[0][1] = EPSILON;
                    let model = get_hyperelastic_constitutive_model();
                    let first_piola_kirchoff_stress = model.calculate_first_piola_kirchoff_stress(&deformation_gradient);
                    assert!((EPSILON*model.get_shear_modulus()/first_piola_kirchoff_stress[0][1] - 1.0).abs() < EPSILON)
                }
                #[test]
                fn zero()
                {
                    get_hyperelastic_constitutive_model().calculate_first_piola_kirchoff_stress(&DeformationGradient::identity()).iter()
                    .for_each(|first_piola_kirchoff_stress_i|
                        first_piola_kirchoff_stress_i.iter()
                        .for_each(|first_piola_kirchoff_stress_ij|
                            assert_eq!(first_piola_kirchoff_stress_ij, &0.0)
                        )
                    )
                }
            }
        }
    }
}
pub(crate) use test_hyperelastic_constitutive_model;