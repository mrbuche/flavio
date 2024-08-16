use super::*;
use super::super::test::*;

use_elastic_hyperviscous_macros!();

test_solid_hyperviscoelastic_constitutive_model!(
    SaintVenantKirchoff, SAINTVENANTKIRCHOFFPARAMETERS,
    SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS)
);

mod consistency
{
    use crate::
    {
        ABS_TOL,
        constitutive::solid::
        {
            elastic::Elastic,
            hyperelastic::
            {
                Hyperelastic,
                SaintVenantKirchoff as HyperelasticSaintVenantKirchoff,
                test::SAINTVENANTKIRCHOFFPARAMETERS as HYPERELASTICSAINTVENANTKIRCHOFFPARAMETERS
            }
        }
    };
    use super::*;
    #[test]
    fn helmholtz_free_energy_density()
    {
        let model = SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS);
        let hyperelastic_model = HyperelasticSaintVenantKirchoff::new(HYPERELASTICSAINTVENANTKIRCHOFFPARAMETERS);
        assert!((model.calculate_helmholtz_free_energy_density(
                &get_deformation_gradient()
            ).unwrap() - hyperelastic_model.calculate_helmholtz_free_energy_density(
                &get_deformation_gradient()
            ).unwrap()).abs() < ABS_TOL
        )
    }
    #[test]
    fn cauchy_stress()
    {
        let model = SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS);
        let hyperelastic_model = HyperelasticSaintVenantKirchoff::new(HYPERELASTICSAINTVENANTKIRCHOFFPARAMETERS);
        model.calculate_cauchy_stress(&get_deformation_gradient(), &DeformationGradientRate::zero()).iter()
        .zip(hyperelastic_model.calculate_cauchy_stress(&get_deformation_gradient()).iter())
        .for_each(|(cauchy_stress_i, elastic_cauchy_stress_i)|
            cauchy_stress_i.iter()
            .zip(elastic_cauchy_stress_i.iter())
            .for_each(|(cauchy_stress_ij, elastic_cauchy_stress_ij)|
                assert!((cauchy_stress_ij - elastic_cauchy_stress_ij).abs() < ABS_TOL)
            )
        )
    }
}
