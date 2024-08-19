use super::super::test::*;
use super::*;

use_thermoelastic_macros!();

test_solid_thermohyperelastic_constitutive_model!(
    SaintVenantKirchoff,
    SAINTVENANTKIRCHOFFPARAMETERS,
    SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS)
);

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
                .calculate_helmholtz_free_energy_density(
                    &get_deformation_gradient(),
                    model.get_reference_temperature()
                )
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
                model.get_reference_temperature(),
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
    #[test]
    fn cauchy_tangent_stiffness() {
        let model = SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS);
        let hyperelastic_model =
            HyperelasticSaintVenantKirchoff::new(HYPERELASTICSAINTVENANTKIRCHOFFPARAMETERS);
        model.calculate_cauchy_tangent_stiffness(&get_deformation_gradient(), model.get_reference_temperature()).expect("the unexpected").iter()
        .zip(hyperelastic_model.calculate_cauchy_tangent_stiffness(&get_deformation_gradient()).expect("the unexpected").iter())
        .for_each(|(cauchy_tangent_stiffness_i, elastic_cauchy_tangent_stiffness_i)|
            cauchy_tangent_stiffness_i.iter()
            .zip(elastic_cauchy_tangent_stiffness_i.iter())
            .for_each(|(cauchy_tangent_stiffness_ij, elastic_cauchy_tangent_stiffness_ij)|
                cauchy_tangent_stiffness_ij.iter()
                .zip(elastic_cauchy_tangent_stiffness_ij.iter())
                .for_each(|(cauchy_tangent_stiffness_ijk, elastic_cauchy_tangent_stiffness_ijk)|
                    cauchy_tangent_stiffness_ijk.iter()
                    .zip(elastic_cauchy_tangent_stiffness_ijk.iter())
                    .for_each(|(cauchy_tangent_stiffness_ijkl, elastic_cauchy_tangent_stiffness_ijkl)|
                        assert!((cauchy_tangent_stiffness_ijkl - elastic_cauchy_tangent_stiffness_ijkl).abs() < ABS_TOL)
                    )
                )
            )
        )
    }
}
