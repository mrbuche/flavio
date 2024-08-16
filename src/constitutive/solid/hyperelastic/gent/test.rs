use super::super::test::*;
use super::*;

use_elastic_macros!();

test_solid_hyperelastic_constitutive_model!(Gent, GENTPARAMETERS, Gent::new(GENTPARAMETERS));

#[test]
fn get_extensibility() {
    assert_eq!(
        &GENTPARAMETERS[2],
        Gent::new(GENTPARAMETERS).get_extensibility()
    )
}

mod panic {
    use super::*;
    #[test]
    #[should_panic]
    fn calculate_cauchy_stress() {
        Gent::new(GENTPARAMETERS).calculate_cauchy_stress(&DeformationGradient::new([
            [GENTPARAMETERS[2], 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]));
    }
    #[test]
    #[should_panic]
    fn calculate_cauchy_tangent_stiffness() {
        Gent::new(GENTPARAMETERS).calculate_cauchy_tangent_stiffness(&DeformationGradient::new([
            [GENTPARAMETERS[2], 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]));
    }
}
