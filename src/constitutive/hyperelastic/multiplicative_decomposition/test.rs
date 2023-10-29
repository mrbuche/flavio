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
    mechanics::test::get_deformation_gradient
};
use super::*;

mod dual
{
    use super::*;
    test_hyperelastic_constitutive_model_constructed!(
        CompositeHyperelasticConstitutiveModelMultiplicativeDecomposition::construct(
            NeoHookeanModel::new(NEOHOOKEANPARAMETERS),
            NeoHookeanModel::new(NEOHOOKEANPARAMETERS)
        )
    );
}

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