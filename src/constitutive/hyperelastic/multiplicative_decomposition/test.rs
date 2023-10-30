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

mod objective
{
    use super::*;
    fn get_composite_constitutive_model<'a>() -> CompositeHyperelasticConstitutiveModelMultiplicativeDecomposition<GentModel<'a>, YeohModel<'a>>
    {
        CompositeHyperelasticConstitutiveModelMultiplicativeDecomposition::construct(
            GentModel::new(GENTPARAMETERS),
            YeohModel::new(YEOHPARAMETERS)
        )
    }
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
                &model.calculate_objective(&get_deformation_gradient(), &get_deformation_gradient_2()),
                &model.calculate_objective(&get_deformation_gradient_rotated(), &get_deformation_gradient_2_rotated())
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
            assert_eq!(get_composite_constitutive_model().calculate_objective(&DeformationGradient::identity(), &DeformationGradient::identity().into()), 0.0)
        }
    }
}