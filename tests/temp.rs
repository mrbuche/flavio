#[cfg(feature = "constitutive")]
mod public {
    mod hyperviscoelastic {
        use flavio::{
            constitutive::{
                hybrid::{Additive, Hybrid, Multiplicative},
                solid::{
                    hyperelastic::{
                        ArrudaBoyce, Gent, Hyperelastic,
                        SaintVenantKirchoff as SaintVenantKirchoff1,
                    },
                    hyperviscoelastic::{
                        Hyperviscoelastic, SaintVenantKirchoff as SaintVenantKirchoff2,
                    },
                    thermohyperelastic::{
                        SaintVenantKirchoff as SaintVenantKirchoff3, Thermohyperelastic,
                    },
                },
                Constitutive,
            },
            mechanics::DeformationGradient,
        };
        #[test]
        #[should_panic]
        fn additive() {
            let model = Additive::construct(
                ArrudaBoyce::new(&[13.0, 3.0, 8.0]),
                Gent::new(&[13.0, 3.0, 23.0]),
            );
            let mut deformation_gradient = DeformationGradient::zero();
            deformation_gradient[0][0] = 8.0;
            (1..3)
                .for_each(|i| deformation_gradient[i][i] = 1.0 / deformation_gradient[0][0].sqrt());
            let result = model.calculate_helmholtz_free_energy_density(&deformation_gradient);
            match result {
                Ok(helmholtz_free_energy_density) => println!("{}", helmholtz_free_energy_density),
                Err(why) => panic!("{}", why),
            }
        }
        #[test]
        #[should_panic]
        fn multiplicative() {
            let model = Multiplicative::construct(
                ArrudaBoyce::new(&[13.0, 3.0, 8.0]),
                Gent::new(&[13.0, 3.0, 23.0]),
            );
            let mut deformation_gradient = DeformationGradient::zero();
            deformation_gradient[0][0] = 8.0;
            (1..3)
                .for_each(|i| deformation_gradient[i][i] = 1.0 / deformation_gradient[0][0].sqrt());
            let result = model.calculate_helmholtz_free_energy_density(&deformation_gradient);
            match result {
                Ok(helmholtz_free_energy_density) => println!("{}", helmholtz_free_energy_density),
                Err(why) => panic!("{}", why),
            }
        }
        #[test]
        #[should_panic]
        fn arruda_boyce() {
            let model = ArrudaBoyce::new(&[13.0, 3.0, 8.0]);
            let mut deformation_gradient = DeformationGradient::zero();
            deformation_gradient[0][0] = 8.0;
            (1..3)
                .for_each(|i| deformation_gradient[i][i] = 1.0 / deformation_gradient[0][0].sqrt());
            let result = model.calculate_helmholtz_free_energy_density(&deformation_gradient);
            match result {
                Ok(helmholtz_free_energy_density) => println!("{}", helmholtz_free_energy_density),
                Err(why) => panic!("{}", why),
            }
        }
        #[test]
        #[should_panic]
        fn gent() {
            let model = Gent::new(&[13.0, 3.0, 23.0]);
            let mut deformation_gradient = DeformationGradient::zero();
            deformation_gradient[0][0] = 8.0;
            (1..3)
                .for_each(|i| deformation_gradient[i][i] = 1.0 / deformation_gradient[0][0].sqrt());
            let result = model.calculate_helmholtz_free_energy_density(&deformation_gradient);
            match result {
                Ok(helmholtz_free_energy_density) => println!("{}", helmholtz_free_energy_density),
                Err(why) => panic!("{}", why),
            }
        }
        #[test]
        #[should_panic]
        fn saint_venant_kirchoff_1() {
            let model = SaintVenantKirchoff1::new(&[1.0, 1.0]);
            let mut deformation_gradient = DeformationGradient::zero();
            (0..3).for_each(|i| deformation_gradient[i][i] = -1.23);
            let result = model.calculate_helmholtz_free_energy_density(&deformation_gradient);
            match result {
                Ok(helmholtz_free_energy_density) => println!("{}", helmholtz_free_energy_density),
                Err(why) => panic!("{}", why),
            }
        }
        #[test]
        #[should_panic]
        fn saint_venant_kirchoff_2() {
            let model = SaintVenantKirchoff2::new(&[1.0, 1.0, 1.0, 1.0]);
            let mut deformation_gradient = DeformationGradient::zero();
            (0..3).for_each(|i| deformation_gradient[i][i] = -1.23);
            let result = model.calculate_helmholtz_free_energy_density(&deformation_gradient);
            match result {
                Ok(helmholtz_free_energy_density) => println!("{}", helmholtz_free_energy_density),
                Err(why) => panic!("{}", why),
            }
        }
        #[test]
        #[should_panic]
        fn saint_venant_kirchoff_3() {
            let model = SaintVenantKirchoff3::new(&[1.0, 1.0, 1.0, 1.0]);
            let mut deformation_gradient = DeformationGradient::zero();
            (0..3).for_each(|i| deformation_gradient[i][i] = -1.23);
            let result =
                model.calculate_helmholtz_free_energy_density(&deformation_gradient, &123.0);
            match result {
                Ok(helmholtz_free_energy_density) => println!("{}", helmholtz_free_energy_density),
                Err(why) => panic!("{}", why),
            }
        }
    }
}
