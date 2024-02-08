use super::
{
    *, super::test::*
};
use crate::
{
    ABS_TOL,
    constitutive::solid::elastic::
    {
        Elastic,
        AlmansiHamel as ElasticAlmansiHamel,
        test::
        {
            ALMANSIHAMELPARAMETERS as ELASTICALMANSIHAMELPARAMETERS
        }
    },
    mechanics::test::get_deformation_gradient
};

crate::constitutive::solid::elastic::test::test_solid_constitutive_model!(
    AlmansiHamel,
    ALMANSIHAMELPARAMETERS,
    AlmansiHamel::new(ALMANSIHAMELPARAMETERS)
);

crate::constitutive::solid::elastic::test::test_solid_elastic_constitutive_model!(
    AlmansiHamel::new(ALMANSIHAMELPARAMETERS)
);



test_thermoelastic_constitutive_model!(
    AlmansiHamel,
    ALMANSIHAMELPARAMETERS,
    AlmansiHamel::new(ALMANSIHAMELPARAMETERS)
);

test_thermoelastic_only_constitutive_model_constructed!(
    AlmansiHamel::new(ALMANSIHAMELPARAMETERS)
);

#[test]
fn cauchy_stress()
{
    let model = AlmansiHamel::new(ALMANSIHAMELPARAMETERS);
    let elastic_model = ElasticAlmansiHamel::new(ELASTICALMANSIHAMELPARAMETERS);
    model.calculate_cauchy_stress(&get_deformation_gradient(), &model.get_reference_temperature()).iter()
    .zip(elastic_model.calculate_cauchy_stress(&get_deformation_gradient()).iter())
    .for_each(|(cauchy_stress_i, elastic_cauchy_stress_i)|
        cauchy_stress_i.iter()
        .zip(elastic_cauchy_stress_i.iter())
        .for_each(|(cauchy_stress_ij, elastic_cauchy_stress_ij)|
            assert!((cauchy_stress_ij - elastic_cauchy_stress_ij).abs() < ABS_TOL)
        )
    )
}

#[test]
fn cauchy_tangent_stiffness()
{
    let model = AlmansiHamel::new(ALMANSIHAMELPARAMETERS);
    let elastic_model = ElasticAlmansiHamel::new(ELASTICALMANSIHAMELPARAMETERS);
    model.calculate_cauchy_tangent_stiffness(&get_deformation_gradient(), &model.get_reference_temperature()).iter()
    .zip(elastic_model.calculate_cauchy_tangent_stiffness(&get_deformation_gradient()).iter())
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