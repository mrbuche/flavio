#[cfg(feature = "constitutive")]
mod public
{
    use flavio::constitutive::
    {
        GentModel,
        NeoHookeanModel,
        YeohModel
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
    #[test]
    fn yeoh_model()
    {
        let _: YeohModel;
    }
}