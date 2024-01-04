#[cfg(feature = "constitutive")]
mod public
{
    use flavio::constitutive::
    {
        elastic::AlmansiHamelModel,
        hyperelastic::
        {
            ArrudaBoyceModel,
            GentModel,
            MooneyRivlinModel,
            NeoHookeanModel,
            YeohModel
        }
    };
    #[test]
    fn almansi_hamel()
    {
        let _: AlmansiHamelModel;
    }
    #[test]
    fn arruda_boyce_model()
    {
        let _: ArrudaBoyceModel;
    }
    #[test]
    fn gent_model()
    {
        let _: GentModel;
    }
    #[test]
    fn mooney_rivlin_model()
    {
        let _: MooneyRivlinModel;
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