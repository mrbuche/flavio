#[cfg(feature = "constitutive")]
mod public
{
    use flavio::constitutive::
    {
        NeoHookeanModel
    };
    #[test]
    fn tensor_rank_0()
    {
        let _: NeoHookeanModel;
    }
}