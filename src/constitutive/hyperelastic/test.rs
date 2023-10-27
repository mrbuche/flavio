macro_rules! test_hyperelastic_constitutive_model
{
    ($hyperelastic_constitutive_model: ident, $constitutive_model_parameters: expr) =>
    {
        fn get_hyperelastic_constitutive_model<'a>() -> $hyperelastic_constitutive_model<'a>
        {
            $hyperelastic_constitutive_model::new($constitutive_model_parameters)
        }
        #[test]
        fn todo()
        {
            todo!("Set up automatic testing, organize properly, include things like objectivity, and test one thing at a time.")
        }
    }
}
pub(crate) use test_hyperelastic_constitutive_model;