use crate::
{
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

#[test]
fn todo()
{
    todo!(
        "Test objective, residual, etc. over finite difference, objectivity, etc. before moving on to full tests below.
        
        Can you mock a model with the needed functions to test using the macro?"
    )
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

mod objective
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
            todo!()
        }
        #[test]
        fn zero()
        {
            assert_eq!(
                get_composite_constitutive_model().calculate_objective(
                    &DeformationGradient::identity(),
                    &DeformationGradient::identity().into()
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
            todo!()
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
            todo!()
        }
        #[test]
        fn zero()
        {
            get_composite_constitutive_model().calculate_residual(
                &DeformationGradient::identity(),
                &DeformationGradient::identity().into()
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
                &DeformationGradient::identity().into()
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