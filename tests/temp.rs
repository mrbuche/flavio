#[cfg(feature = "constitutive")]
mod public
{
    mod hyperviscoelastic
    {
        use flavio::
        {
            constitutive::
            {
                solid::hyperviscoelastic::*,
                Constitutive
            },
            mechanics::DeformationGradient
        };
        #[test]
        fn saint_venant_kirchoff()
        {
            let model = SaintVenantKirchoff::new(&[1.0, 1.0, 1.0, 1.0]);
            let mut deformation_gradient = DeformationGradient::zero();
            (0..3).for_each(|i| deformation_gradient[i][i] = -1.23);
            let result = model.calculate_helmholtz_free_energy_density(&deformation_gradient);
            match result {
                Ok(helmholtz_free_energy_density) => println!("{}", helmholtz_free_energy_density),
                Err(why) => panic!("{}", why)
            }
        }
    }
}
