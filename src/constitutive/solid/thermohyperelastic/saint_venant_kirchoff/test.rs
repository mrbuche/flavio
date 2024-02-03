use super::
{
    SaintVenantKirchoff,
    Thermohyperelastic,
    super::test::
    {
        SAINTVENANTKIRCHOFFPARAMETERS,
        test_thermohyperelastic_constitutive_model
    }
};
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
    },
    mechanics::test::get_deformation_gradient
};

test_thermohyperelastic_constitutive_model!(
    SaintVenantKirchoff,
    SAINTVENANTKIRCHOFFPARAMETERS,
    SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS)
);

#[test]
fn helmholtz_free_energy_density()
{
    let model = SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS);
    let hyperelastic_model = HyperelasticSaintVenantKirchoff::new(HYPERELASTICSAINTVENANTKIRCHOFFPARAMETERS);
    assert!(
        (model.calculate_helmholtz_free_energy_density(&get_deformation_gradient(), &model.get_reference_temperature())
        - hyperelastic_model.calculate_helmholtz_free_energy_density(&get_deformation_gradient())).abs() < ABS_TOL
    )
}

#[test]
fn cauchy_stress()
{
    let model = SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS);
    let hyperelastic_model = HyperelasticSaintVenantKirchoff::new(HYPERELASTICSAINTVENANTKIRCHOFFPARAMETERS);
    model.calculate_cauchy_stress(&get_deformation_gradient(), &model.get_reference_temperature()).iter()
    .zip(hyperelastic_model.calculate_cauchy_stress(&get_deformation_gradient()).iter())
    .for_each(|(cauchy_stress_i, hyperelastic_cauchy_stress_i)|
        cauchy_stress_i.iter()
        .zip(hyperelastic_cauchy_stress_i.iter())
        .for_each(|(cauchy_stress_ij, hyperelastic_cauchy_stress_ij)|
            assert!((cauchy_stress_ij - hyperelastic_cauchy_stress_ij).abs() < ABS_TOL)
        )
    )
}

#[test]
fn cauchy_tangent_stiffness()
{
    let model = SaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS);
    let hyperelastic_model = HyperelasticSaintVenantKirchoff::new(SAINTVENANTKIRCHOFFPARAMETERS);
    model.calculate_cauchy_tangent_stiffness(&get_deformation_gradient(), &model.get_reference_temperature()).iter()
    .zip(hyperelastic_model.calculate_cauchy_tangent_stiffness(&get_deformation_gradient()).iter())
    .for_each(|(cauchy_tangent_stiffness_i, hyperelastic_cauchy_tangent_stiffness_i)|
        cauchy_tangent_stiffness_i.iter()
        .zip(hyperelastic_cauchy_tangent_stiffness_i.iter())
        .for_each(|(cauchy_tangent_stiffness_ij, hyperelastic_cauchy_tangent_stiffness_ij)|
            cauchy_tangent_stiffness_ij.iter()
            .zip(hyperelastic_cauchy_tangent_stiffness_ij.iter())
            .for_each(|(cauchy_tangent_stiffness_ijk, hyperelastic_cauchy_tangent_stiffness_ijk)|
                cauchy_tangent_stiffness_ijk.iter()
                .zip(hyperelastic_cauchy_tangent_stiffness_ijk.iter())
                .for_each(|(cauchy_tangent_stiffness_ijkl, hyperelastic_cauchy_tangent_stiffness_ijkl)|
                    assert!((cauchy_tangent_stiffness_ijkl - hyperelastic_cauchy_tangent_stiffness_ijkl).abs() < ABS_TOL)
                )
            )
        )
    )
}