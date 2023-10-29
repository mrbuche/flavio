use crate::constitutive::
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
};
use super::*;

mod dual
{
    use super::*;
    test_hyperelastic_constitutive_model_constructed!(
        CompositeHyperelasticConstitutiveModelEqualDeformation::construct(
            NeoHookeanModel::new(NEOHOOKEANPARAMETERS),
            NeoHookeanModel::new(NEOHOOKEANPARAMETERS)
        )
    );
}

mod mized
{
    use super::*;
    test_hyperelastic_constitutive_model_constructed!(
        CompositeHyperelasticConstitutiveModelEqualDeformation::construct(
            GentModel::new(GENTPARAMETERS),
            YeohModel::new(YEOHPARAMETERS)
        )
    );
}

mod nested
{
    use super::*;
    test_hyperelastic_constitutive_model_constructed!(
        CompositeHyperelasticConstitutiveModelEqualDeformation::construct(
            CompositeHyperelasticConstitutiveModelEqualDeformation::construct(
                GentModel::new(GENTPARAMETERS),
                YeohModel::new(YEOHPARAMETERS)
            ),
            CompositeHyperelasticConstitutiveModelEqualDeformation::construct(
                MooneyRivlinModel::new(MOONEYRIVLINPARAMETERS),
                NeoHookeanModel::new(NEOHOOKEANPARAMETERS)
            )
        )
    );
}

#[test]
fn still_need_to_test_things_like_equal_stress_for_same_model()
{
    todo!()
}
