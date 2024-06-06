use super::*;
use super::super::test::*;

use_viscoelastic_macros!();

test_solid_elastic_hyperviscous_constitutive_model!(
    AlmansiHamel, ALMANSIHAMELPARAMETERS,
    AlmansiHamel::new(ALMANSIHAMELPARAMETERS)
);

mod consistency
{
    use crate::
    {
        ABS_TOL,
        constitutive::solid::
        {
            elastic::
            {
                Elastic,
                AlmansiHamel as ElasticAlmansiHamel,
                test::ALMANSIHAMELPARAMETERS as ELASTICALMANSIHAMELPARAMETERS
            }
        }
    };
    use super::*;
    #[test]
    fn cauchy_stress()
    {
        let model = AlmansiHamel::new(ALMANSIHAMELPARAMETERS);
        let hyperelastic_model = ElasticAlmansiHamel::new(ELASTICALMANSIHAMELPARAMETERS);
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
