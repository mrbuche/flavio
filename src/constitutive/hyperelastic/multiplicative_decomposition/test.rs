use crate::
{
    EPSILON,
    constitutive::
    {
        GentModel,
        MooneyRivlinModel,
        NeoHookeanModel,
        YeohModel,
        hyperelastic::test::test_hyperelastic_constitutive_model_constructed,
        test::
        {
            GENTPARAMETERS,
            MOONEYRIVLINPARAMETERS,
            NEOHOOKEANPARAMETERS,
            YEOHPARAMETERS
        }
    },
    mechanics::test::
    {
        get_deformation_gradient,
        get_deformation_gradient_2,
        get_deformation_gradient_rotated,
        get_deformation_gradient_2_rotated,
        get_rotation_current_configuration,
        get_rotation_intermediate_configuration,
        get_rotation_reference_configuration
    },
    test::assert_eq_within_tols
};
use super::*;

#[test]
fn dont_forget_to_test_mandel_stress_too()
{
    todo!()
}

#[test]
fn also_consider_the_combo_methods_returning_tuples()
{
    todo!()
}

// mod dual
// {
//     use super::*;
//     test_hyperelastic_constitutive_model_constructed!(
//         CompositeHyperelasticConstitutiveModelMultiplicativeDecomposition::construct(
//             NeoHookeanModel::new(NEOHOOKEANPARAMETERS),
//             NeoHookeanModel::new(NEOHOOKEANPARAMETERS)
//         )
//     );
// }

fn get_composite_constitutive_model<'a>() -> CompositeHyperelasticConstitutiveModelMultiplicativeDecomposition<GentModel<'a>, YeohModel<'a>>
{
    CompositeHyperelasticConstitutiveModelMultiplicativeDecomposition::construct(
        GentModel::new(GENTPARAMETERS),
        YeohModel::new(YEOHPARAMETERS)
    )
}

fn calculate_residual_from_finite_difference_of_objective(is_deformed: bool) -> FirstPiolaKirchoffStress2
{
    let model = get_composite_constitutive_model();
    let mut residual = FirstPiolaKirchoffStress2::zero();
    let deformation_gradient =
        if is_deformed
        {
            get_deformation_gradient()
        }
        else
        {
            DeformationGradient::identity()
        };
    for i in 0..3
    {
        for j in 0..3
        {
            let mut deformation_gradient_2_plus = 
                if is_deformed
                {
                    get_deformation_gradient_2()
                }
                else
                {
                    DeformationGradient2::identity()
                };
            deformation_gradient_2_plus[i][j] += 0.5*EPSILON;
            let objective_plus = model.calculate_objective(
                &deformation_gradient,
                &deformation_gradient_2_plus
            );
            let mut deformation_gradient_2_minus = 
                if is_deformed
                {
                    get_deformation_gradient_2()
                }
                else
                {
                    DeformationGradient2::identity()
                };
            deformation_gradient_2_minus[i][j] -= 0.5*EPSILON;
            let objective_minus = model.calculate_objective(
                &deformation_gradient,
                &deformation_gradient_2_minus
            );
            residual[i][j] = (objective_plus - objective_minus)/EPSILON;
        }
    }
    residual
}

fn calculate_first_piola_kirchoff_tangent_stiffness_from_finite_difference_of_first_piola_kirchoff_stress(is_deformed: bool) -> FirstPiolaKirchoffTangentStiffness2
{
    let model = get_composite_constitutive_model();
    let mut residual_tangent = FirstPiolaKirchoffTangentStiffness2::zero();
    let deformation_gradient =
        if is_deformed
        {
            get_deformation_gradient()
        }
        else
        {
            DeformationGradient::identity()
        };
    for k in 0..3
    {
        for l in 0..3
        {
            let mut deformation_gradient_2_plus = 
                if is_deformed
                {
                    get_deformation_gradient_2()
                }
                else
                {
                    DeformationGradient2::identity()
                };
            deformation_gradient_2_plus[k][l] += 0.5*EPSILON;
            let residual_plus = model.calculate_residual(
                &deformation_gradient,
                &deformation_gradient_2_plus
            );
            let mut deformation_gradient_2_minus = 
                if is_deformed
                {
                    get_deformation_gradient_2()
                }
                else
                {
                    DeformationGradient2::identity()
                };
            deformation_gradient_2_minus[k][l] -= 0.5*EPSILON;
            let residual_minus = model.calculate_residual(
                &deformation_gradient,
                &deformation_gradient_2_minus
            );
            for i in 0..3
            {
                for j in 0..3
                {
                    residual_tangent[i][j][k][l] = (residual_plus[i][j] - residual_minus[i][j])/EPSILON;
                }
            }
        }
    }
    residual_tangent
}

mod objective
{
    use super::*;
    mod deformed
    {
        use super::*;
        #[test]
        fn finite_difference()
        {
            get_composite_constitutive_model().calculate_residual(
                &get_deformation_gradient(),
                &get_deformation_gradient_2()
            ).iter().zip(
                calculate_residual_from_finite_difference_of_objective(true).iter()
            ).for_each(|(residual_i, fd_residual_i)|
                residual_i.iter()
                .zip(fd_residual_i.iter())
                .for_each(|(residual_ij, fd_residual_ij)|
                    assert!((residual_ij/fd_residual_ij - 1.0).abs() < EPSILON)
                )
            )
        }
        #[test]
        fn objectivity()
        {
            let model = get_composite_constitutive_model();
            assert_eq_within_tols(
                &model.calculate_objective(
                    &get_deformation_gradient(),
                    &get_deformation_gradient_2()
                ), &model.calculate_objective(
                    &get_deformation_gradient_rotated(),
                    &get_deformation_gradient_2_rotated()
                )
            )
        }
    }
    mod undeformed
    {
        use super::*;
        #[test]
        fn finite_difference()
        {
            calculate_residual_from_finite_difference_of_objective(false).iter()
            .for_each(|fd_residual_i|
                fd_residual_i.iter()
                .for_each(|fd_residual_ij|
                    assert!(fd_residual_ij.abs() < EPSILON)
                )
            )
        }
        #[test]
        fn zero()
        {
            assert_eq!(
                get_composite_constitutive_model().calculate_objective(
                    &DeformationGradient::identity(),
                    &DeformationGradient2::identity()
                ), 0.0
            )
        }
    }
}

mod residual
{
    use super::*;
    mod deformed
    {
        use super::*;
        #[test]
        fn finite_difference()
        {
            get_composite_constitutive_model().calculate_residual_tangent(
                &get_deformation_gradient(),
                &get_deformation_gradient_2()
            ).iter().zip(
                calculate_first_piola_kirchoff_tangent_stiffness_from_finite_difference_of_first_piola_kirchoff_stress(true).iter()
            ).for_each(|(residual_tangent_i, fd_residual_tangent_i)|
                residual_tangent_i.iter()
                .zip(fd_residual_tangent_i.iter())
                .for_each(|(residual_tangent_ij, fd_residual_tangent_ij)|
                    residual_tangent_ij.iter()
                    .zip(fd_residual_tangent_ij.iter())
                    .for_each(|(residual_tangent_ijk, fd_residual_tangent_ijk)|
                        residual_tangent_ijk.iter()
                        .zip(fd_residual_tangent_ijk.iter())
                        .for_each(|(residual_tangent_ijkl, fd_residual_tangent_ijkl)|
                            assert!((residual_tangent_ijkl/fd_residual_tangent_ijkl - 1.0).abs() < EPSILON)
                        )
                    )
                )
            )
        }
        #[test]
        fn objectivity()
        {
            let model = get_composite_constitutive_model();
            model.calculate_residual(
                &get_deformation_gradient(),
                &get_deformation_gradient_2()
            ).iter().zip((
                get_rotation_intermediate_configuration().transpose() *
                model.calculate_residual(
                    &get_deformation_gradient_rotated(),
                    &get_deformation_gradient_2_rotated()
                ) * get_rotation_reference_configuration()
            ).iter())
            .for_each(|(residuali, rotated_residual_i)|
                residuali.iter()
                .zip(rotated_residual_i.iter())
                .for_each(|(residualij, rotated_residual_ij)|
                    assert_eq_within_tols(residualij, rotated_residual_ij)
                )
            )
        }
    }
    mod undeformed
    {
        use super::*;
        #[test]
        fn finite_difference()
        {
            get_composite_constitutive_model().calculate_residual_tangent(
                &DeformationGradient::identity(),
                &DeformationGradient2::identity()
            ).iter().zip(
                calculate_first_piola_kirchoff_tangent_stiffness_from_finite_difference_of_first_piola_kirchoff_stress(false).iter()
            ).for_each(|(residual_tangent_i, fd_residual_tangent_i)|
                residual_tangent_i.iter()
                .zip(fd_residual_tangent_i.iter())
                .for_each(|(residual_tangent_ij, fd_residual_tangent_ij)|
                    residual_tangent_ij.iter()
                    .zip(fd_residual_tangent_ij.iter())
                    .for_each(|(residual_tangent_ijk, fd_residual_tangent_ijk)|
                        residual_tangent_ijk.iter()
                        .zip(fd_residual_tangent_ijk.iter())
                        .for_each(|(residual_tangent_ijkl, fd_residual_tangent_ijkl)|
                            assert!(
                                (residual_tangent_ijkl/fd_residual_tangent_ijkl - 1.0).abs() < EPSILON ||
                                fd_residual_tangent_ijkl.abs() < EPSILON
                            )
                        )
                    )
                )
            )
        }
        #[test]
        fn zero()
        {
            get_composite_constitutive_model().calculate_residual(
                &DeformationGradient::identity(),
                &DeformationGradient2::identity()
            ).iter().for_each(|residual_i|
                residual_i.iter().for_each(|residual_ij|
                    assert_eq!(residual_ij, &0.0)
                )
            )
        }
    }
}

mod residual_tangent
{
    use super::*;
    mod deformed
    {
        use super::*;
        #[test]
        fn objectivity()
        {
            let model = get_composite_constitutive_model();
            model.calculate_residual_tangent(
                &get_deformation_gradient(),
                &get_deformation_gradient_2()
            ).iter().zip((
                model.calculate_residual_tangent(
                    &get_deformation_gradient_rotated(),
                    &get_deformation_gradient_2_rotated()
                ).contract_all_indices_with_first_indices_of(
                    &get_rotation_intermediate_configuration(),
                    &get_rotation_reference_configuration(),
                    &get_rotation_intermediate_configuration(),
                    &get_rotation_reference_configuration()
                )
            ).iter())
            .for_each(|(residual_tangent_i, rotated_residual_tangent_i)|
                residual_tangent_i.iter()
                .zip(rotated_residual_tangent_i.iter())
                .for_each(|(residual_tangent_ij, rotated_residual_tangent_ij)|
                    residual_tangent_ij.iter()
                    .zip(rotated_residual_tangent_ij.iter())
                    .for_each(|(residual_tangent_ijk, rotated_residual_tangent_ijk)|
                        residual_tangent_ijk.iter()
                        .zip(rotated_residual_tangent_ijk.iter())
                        .for_each(|(residual_tangent_ijkl, rotated_residual_tangent_ijkl)|
                            assert_eq_within_tols(residual_tangent_ijkl, rotated_residual_tangent_ijkl)
                        )
                    )
                )
            )
        }
        #[test]
        fn symmetry()
        {
            let model = get_composite_constitutive_model();
            let residual_tangent = model.calculate_residual_tangent(
                &get_deformation_gradient(),
                &get_deformation_gradient_2()
            );
            for i in 0..3
            {
                for j in 0..3
                {
                    for k in 0..3
                    {
                        for l in 0..3
                        {
                            assert_eq_within_tols(
                                &residual_tangent[i][j][k][l],
                                &residual_tangent[k][l][i][j]
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
        fn symmetry()
        {
            let model = get_composite_constitutive_model();
            let residual_tangent = model.calculate_residual_tangent(
                &DeformationGradient::identity(),
                &DeformationGradient2::identity()
            );
            for i in 0..3
            {
                for j in 0..3
                {
                    for k in 0..3
                    {
                        for l in 0..3
                        {
                            assert_eq_within_tols(
                                &residual_tangent[i][j][k][l],
                                &residual_tangent[k][l][i][j]
                            )
                        }
                    }
                }
            }
        }
    }
}