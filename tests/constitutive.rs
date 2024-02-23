#[cfg(feature = "constitutive")]
mod public
{
    mod elastic
    {
        use flavio::constitutive::solid::elastic::AlmansiHamel;
        #[test]
        fn almansi_hamel()
        {
            let _: AlmansiHamel;
        }
    }
    mod hyperelastic
    {
        use flavio::constitutive::solid::hyperelastic::
        {
            ArrudaBoyce,
            Gent,
            Fung,
            MooneyRivlin,
            NeoHookean,
            SaintVenantKirchoff,
            Yeoh
        };
        #[test]
        fn arruda_boyce()
        {
            let _: ArrudaBoyce;
        }
        #[test]
        fn fung()
        {
            let _: Fung;
        }
        #[test]
        fn gent()
        {
            let _: Gent;
        }
        #[test]
        fn mooney_rivlin()
        {
            let _: MooneyRivlin;
        }
        #[test]
        fn neo_hookean()
        {
            let _: NeoHookean;
        }
        #[test]
        fn saint_venant_kirchoff()
        {
            let _: SaintVenantKirchoff;
        }
        #[test]
        fn yeoh()
        {
            let _: Yeoh;
        }
    }
    mod elastic_hyperviscous
    {
        use flavio::constitutive::solid::elastic_hyperviscous::AlmansiHamel;
        #[test]
        fn almansi_hamel()
        {
            let _: AlmansiHamel;
        }
    }
    mod hyperviscoelastic
    {
        use flavio::constitutive::solid::hyperviscoelastic::SaintVenantKirchoff;
        #[test]
        fn saint_venant_kirchoff()
        {
            let _: SaintVenantKirchoff;
        }
    }
    mod thermoelastic
    {
        use flavio::constitutive::solid::thermoelastic::AlmansiHamel;
        #[test]
        fn almansi_hamel()
        {
            let _: AlmansiHamel;
        }
    }
    mod thermohyperelastic
    {
        use flavio::constitutive::solid::thermohyperelastic::SaintVenantKirchoff;
        #[test]
        fn almansi_hamel()
        {
            let _: SaintVenantKirchoff;
        }
    }
}