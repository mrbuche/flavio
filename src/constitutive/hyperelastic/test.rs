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
            math::TensorRank2Trait,
            mechanics::test::
            {
                get_deformation_gradient,
                get_deformation_gradient_rotated,
                get_rotation_current_configuration,
                get_rotation_reference_configuration
            },
            test::assert_eq_within_tols
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
        mod cauchy_stress
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    todo!()
                }
                #[test]
                fn objectivity()
                {
                    let model = get_hyperelastic_constitutive_model();
                    model.calculate_cauchy_stress(&get_deformation_gradient()).iter()
                    .zip((
                        get_rotation_current_configuration().transpose() *
                        model.calculate_cauchy_stress(&get_deformation_gradient_rotated()) *
                        get_rotation_current_configuration()
                    ).iter())
                    .for_each(|(cauchy_stress_i, rotated_cauchy_stress_i)|
                        cauchy_stress_i.iter()
                        .zip(rotated_cauchy_stress_i.iter())
                        .for_each(|(cauchy_stress_ij, rotated_cauchy_stress_ij)|
                            assert_eq_within_tols(cauchy_stress_ij, rotated_cauchy_stress_ij)
                        )
                    )
                }
                #[test]
                fn symmetry()
                {
                    let model = get_hyperelastic_constitutive_model();
                    let cauchy_stress = model.calculate_cauchy_stress(&get_deformation_gradient());
                    for i in 0..3
                    {
                        for j in 0..3
                        {
                            assert_eq_within_tols(
                                &cauchy_stress[i][j],
                                &cauchy_stress[j][i]
                            )
                        }
                    }
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn bulk_modulus()
                {
                    let model = get_hyperelastic_constitutive_model();
                    let deformation_gradient = DeformationGradient::identity() * (1.0 + EPSILON).powf(1.0/3.0);
                    let cauchy_stress = model.calculate_cauchy_stress(&deformation_gradient);
                    assert!((3.0*EPSILON*model.get_bulk_modulus()/cauchy_stress.trace() - 1.0).abs() < EPSILON)
                }
                #[test]
                fn finite_difference()
                {
                    todo!()
                }
                #[test]
                fn shear_modulus()
                {
                    let model = get_hyperelastic_constitutive_model();
                    let mut deformation_gradient = DeformationGradient::identity();
                    deformation_gradient[0][1] = EPSILON;
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
        mod cauchy_tangent_stiffness
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    todo!()
                }
                #[test]
                fn objectivity()
                {
                    todo!()
                }
                #[test]
                fn symmetry()
                {
                    let model = get_hyperelastic_constitutive_model();
                    let cauchy_tangent_stiffness = model.calculate_cauchy_tangent_stiffness(&get_deformation_gradient());
                    for i in 0..3
                    {
                        for j in 0..3
                        {
                            for k in 0..3
                            {
                                for l in 0..3
                                {
                                    assert_eq_within_tols(
                                        &cauchy_tangent_stiffness[i][j][k][l],
                                        &cauchy_tangent_stiffness[j][i][k][l]
                                    )
                                }
                            }
                        }
                    }
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    todo!()
                }
                #[test]
                fn symmetry()
                {
                    let model = get_hyperelastic_constitutive_model();
                    let cauchy_tangent_stiffness = model.calculate_cauchy_tangent_stiffness(&DeformationGradient::identity());
                    for i in 0..3
                    {
                        for j in 0..3
                        {
                            for k in 0..3
                            {
                                for l in 0..3
                                {
                                    assert_eq!(cauchy_tangent_stiffness[i][j][k][l], cauchy_tangent_stiffness[j][i][k][l])
                                }
                            }
                        }
                    }
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
                fn finite_difference()
                {
                    todo!()
                }
                #[test]
                fn objectivity()
                {
                    let model = get_hyperelastic_constitutive_model();
                    assert_eq_within_tols(
                        &model.calculate_helmholtz_free_energy_density(&get_deformation_gradient()),
                        &model.calculate_helmholtz_free_energy_density(&get_deformation_gradient_rotated())
                    )
                }
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
                fn finite_difference()
                {
                    todo!()
                }
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
                #[test]
                fn finite_difference()
                {
                    todo!()
                }
                #[test]
                fn objectivity()
                {
                    let model = get_hyperelastic_constitutive_model();
                    model.calculate_first_piola_kirchoff_stress(&get_deformation_gradient()).iter()
                    .zip((
                        get_rotation_current_configuration().transpose() *
                        model.calculate_first_piola_kirchoff_stress(&get_deformation_gradient_rotated()) *
                        get_rotation_reference_configuration()
                    ).iter())
                    .for_each(|(first_piola_kirchoff_stress_i, rotated_first_piola_kirchoff_stress_i)|
                        first_piola_kirchoff_stress_i.iter()
                        .zip(rotated_first_piola_kirchoff_stress_i.iter())
                        .for_each(|(first_piola_kirchoff_stress_ij, rotated_first_piola_kirchoff_stress_ij)|
                            assert_eq_within_tols(first_piola_kirchoff_stress_ij, rotated_first_piola_kirchoff_stress_ij)
                        )
                    )
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn bulk_modulus()
                {
                    let model = get_hyperelastic_constitutive_model();
                    let deformation_gradient = DeformationGradient::identity() * (1.0 + EPSILON).powf(1.0/3.0);
                    let first_piola_kirchoff_stress = model.calculate_first_piola_kirchoff_stress(&deformation_gradient);
                    assert!((3.0*EPSILON*model.get_bulk_modulus()/first_piola_kirchoff_stress.trace() - 1.0).abs() < EPSILON)
                }
                #[test]
                fn finite_difference()
                {
                    todo!()
                }
                #[test]
                fn shear_modulus()
                {
                    let model = get_hyperelastic_constitutive_model();
                    let mut deformation_gradient = DeformationGradient::identity();
                    deformation_gradient[0][1] = EPSILON;
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
        mod first_piola_kirchoff_tangent_stiffness
        {
            use super::*;
            mod deformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    todo!()
                }
                #[test]
                fn objectivity()
                {
                    todo!()
                }
                #[test]
                fn symmetry()
                {
                    let model = get_hyperelastic_constitutive_model();
                    let first_piola_kirchoff_tangent_stiffness = model.calculate_first_piola_kirchoff_tangent_stiffness(&get_deformation_gradient());
                    for i in 0..3
                    {
                        for j in 0..3
                        {
                            for k in 0..3
                            {
                                for l in 0..3
                                {
                                    assert_eq_within_tols(
                                        &first_piola_kirchoff_tangent_stiffness[i][j][k][l],
                                        &first_piola_kirchoff_tangent_stiffness[k][l][i][j]
                                    )
                                }
                            }
                        }
                    }
                }
            }
            mod undeformed
            {
                use super::*;
                #[test]
                fn finite_difference()
                {
                    todo!()
                }
                #[test]
                fn symmetry()
                {
                    let model = get_hyperelastic_constitutive_model();
                    let first_piola_kirchoff_tangent_stiffness = model.calculate_first_piola_kirchoff_tangent_stiffness(&DeformationGradient::identity());
                    for i in 0..3
                    {
                        for j in 0..3
                        {
                            for k in 0..3
                            {
                                for l in 0..3
                                {
                                    assert_eq_within_tols(
                                        &first_piola_kirchoff_tangent_stiffness[i][j][k][l],
                                        &first_piola_kirchoff_tangent_stiffness[k][l][i][j]
                                    )
                                }
                            }
                        }
                    }
                }
            }

        }
    }
}
pub(crate) use test_hyperelastic_constitutive_model;