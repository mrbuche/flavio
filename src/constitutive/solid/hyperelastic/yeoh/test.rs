use super::
{
    YeohModel,
    super::test::
    {
        YEOHPARAMETERS,
        test_hyperelastic_constitutive_model
    }
};

test_hyperelastic_constitutive_model!(
    YeohModel,
    YEOHPARAMETERS,
    YeohModel::new(YEOHPARAMETERS)
);

#[test]
fn get_moduli()
{
    YeohModel::new(YEOHPARAMETERS).get_moduli().iter()
    .zip(YEOHPARAMETERS[1..].iter())
    .for_each(|(modulus_i, parameter_i)|
        assert_eq!(modulus_i, parameter_i)
    )
}

#[test]
fn get_extra_moduli()
{
    YeohModel::new(YEOHPARAMETERS).get_extra_moduli().iter()
    .zip(YEOHPARAMETERS[2..].iter())
    .for_each(|(modulus_i, parameter_i)|
        assert_eq!(modulus_i, parameter_i)
    )
}

