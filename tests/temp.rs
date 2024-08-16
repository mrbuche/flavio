#[cfg(feature = "constitutive")]
mod public
{
    mod hyperviscoelastic
    {
        use flavio::
        {
            constitutive::
            {
                solid::
                {
                    hyperviscoelastic::
                    {
                        SaintVenantKirchoff as SaintVenantKirchoff1,
                        Hyperviscoelastic
                    },
                    thermohyperelastic::
                    {
                        SaintVenantKirchoff as SaintVenantKirchoff2,
                        Thermohyperelastic
                    },
                },
                Constitutive
            },
            mechanics::DeformationGradient
        };
        #[test]
        fn saint_venant_kirchoff_1()
        {
            let model = SaintVenantKirchoff1::new(&[1.0, 1.0, 1.0, 1.0]);
            let mut deformation_gradient = DeformationGradient::zero();
            (0..3).for_each(|i| deformation_gradient[i][i] = -1.23);
            let result = model.calculate_helmholtz_free_energy_density(&deformation_gradient);
            match result {
                Ok(helmholtz_free_energy_density) => println!("{}", helmholtz_free_energy_density),
                Err(why) => panic!("{}", why)
            }
        }
        #[test]
        fn saint_venant_kirchoff_2()
        {
            let model = SaintVenantKirchoff2::new(&[1.0, 1.0, 1.0, 1.0]);
            let mut deformation_gradient = DeformationGradient::zero();
            (0..3).for_each(|i| deformation_gradient[i][i] = -1.23);
            let result = model.calculate_helmholtz_free_energy_density(&deformation_gradient, &123.0);
            match result {
                Ok(helmholtz_free_energy_density) => println!("{}", helmholtz_free_energy_density),
                Err(why) => panic!("{}", why)
            }
        }
    }
}
