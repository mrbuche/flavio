#[cfg(feature = "constitutive")]
mod public
{
    use flavio::constitutive::
    {
        GentModel,
        NeoHookeanModel
    };
    #[test]
    fn gent_model()
    {
        let _: GentModel;
    }
    #[test]
    fn neo_hookean_model()
    {
        let _: NeoHookeanModel;
    }
}