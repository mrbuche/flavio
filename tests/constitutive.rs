#[cfg(feature = "constitutive")]
mod public
{
    use flavio::constitutive::
    {
        NeoHookeanModel
    };
    #[test]
    fn neo_hookean_model()
    {
        let _: NeoHookeanModel;
    }
}