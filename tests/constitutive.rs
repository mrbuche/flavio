#[cfg(feature = "constitutive")]
mod public
{
    use flavio::constitutive::solid::
    {
        elastic::AlmansiHamel,
        hyperelastic::
        {
            ArrudaBoyce,
            Gent,
            Fung,
            MooneyRivlin,
            NeoHookean,
            SaintVenantKirchoff,
            Yeoh
        }
    };
    #[test]
    fn almansi_hamel()
    {
        let _: AlmansiHamel;
    }
    #[test]
    fn arruda_boyce_model()
    {
        let _: ArrudaBoyce;
    }
    #[test]
    fn fung_model()
    {
        let _: Fung;
    }
    #[test]
    fn gent_model()
    {
        let _: Gent;
    }
    #[test]
    fn mooney_rivlin_model()
    {
        let _: MooneyRivlin;
    }
    #[test]
    fn neo_hookean_model()
    {
        let _: NeoHookean;
    }
    #[test]
    fn saint_venant_kirchoff_model()
    {
        let _: SaintVenantKirchoff;
    }
    #[test]
    fn yeoh_model()
    {
        let _: Yeoh;
    }
}