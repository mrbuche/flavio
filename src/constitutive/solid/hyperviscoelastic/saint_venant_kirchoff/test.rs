use super::super::test::*;
use super::*;

use_elastic_hyperviscous_macros!();

test_solid_hyperviscoelastic_constitutive_model!(
    SaintVenantKirchoff,
    SAINTVENANTKIRCHOFFPARAMETERS,
    SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS)
);

test_solve!(SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS));

mod consistency {
    use super::*;
    use crate::{
        constitutive::solid::{
            elastic::Elastic,
            hyperelastic::{
                test::SAINTVENANTKIRCHOFFPARAMETERS as HYPERELASTICSAINTVENANTKIRCHOFFPARAMETERS,
                Hyperelastic, SaintVenantKirchoff as HyperelasticSaintVenantKirchoff,
            },
        },
        ABS_TOL,
    };
    #[test]
    fn helmholtz_free_energy_density() {
        let model = SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS);
        let hyperelastic_model =
            HyperelasticSaintVenantKirchoff::new(HYPERELASTICSAINTVENANTKIRCHOFFPARAMETERS);
        assert!(
            (model
                .calculate_helmholtz_free_energy_density(&get_deformation_gradient())
                .expect("the unexpected")
                - hyperelastic_model
                    .calculate_helmholtz_free_energy_density(&get_deformation_gradient())
                    .expect("the unexpected"))
            .abs()
                < ABS_TOL
        )
    }
    #[test]
    fn cauchy_stress() {
        let model = SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS);
        let hyperelastic_model =
            HyperelasticSaintVenantKirchoff::new(HYPERELASTICSAINTVENANTKIRCHOFFPARAMETERS);
        model
            .calculate_cauchy_stress(
                &get_deformation_gradient(),
                &DeformationGradientRate::zero(),
            )
            .expect("the unexpected")
            .iter()
            .zip(
                hyperelastic_model
                    .calculate_cauchy_stress(&get_deformation_gradient())
                    .expect("the unexpected")
                    .iter(),
            )
            .for_each(|(cauchy_stress_i, elastic_cauchy_stress_i)| {
                cauchy_stress_i
                    .iter()
                    .zip(elastic_cauchy_stress_i.iter())
                    .for_each(|(cauchy_stress_ij, elastic_cauchy_stress_ij)| {
                        assert!((cauchy_stress_ij - elastic_cauchy_stress_ij).abs() < ABS_TOL)
                    })
            })
    }
}
