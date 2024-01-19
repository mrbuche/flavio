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
            FungModel,
            MooneyRivlinModel,
            NeoHookeanModel,
            SaintVenantKirchoffModel,
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
    fn fung_model()
    {
        let _: FungModel;
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
    fn saint_venant_kirchoff_model()
    {
        let _: SaintVenantKirchoffModel;
    }
    #[test]
    fn yeoh_model()
    {
        let _: YeohModel;
    }
}