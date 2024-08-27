use super::super::test::*;
use super::*;

use_elastic_macros!();

test_solid_hyperelastic_constitutive_model!(Gent, GENTPARAMETERS, Gent::new(GENTPARAMETERS));

test_solve!(Gent::new(GENTPARAMETERS));

#[test]
fn get_extensibility() {
    assert_eq!(
        &GENTPARAMETERS[2],
        Gent::new(GENTPARAMETERS).get_extensibility()
    )
}

mod maximum_extensibility {
    use super::*;
    #[test]
    #[should_panic(expected = "Maximum extensibility reached.")]
    fn calculate_cauchy_stress() {
        Gent::new(GENTPARAMETERS)
            .calculate_cauchy_stress(&DeformationGradient::new([
                [16.0, 0.00, 0.00],
                [0.0, 0.25, 0.00],
                [0.0, 0.00, 0.25],
            ]))
            .expect("the unexpected");
    }
    #[test]
    #[should_panic(expected = "Maximum extensibility reached.")]
    fn calculate_cauchy_tangent_stiffness() {
        Gent::new(GENTPARAMETERS)
            .calculate_cauchy_tangent_stiffness(&DeformationGradient::new([
                [16.0, 0.00, 0.00],
                [0.0, 0.25, 0.00],
                [0.0, 0.00, 0.25],
            ]))
            .expect("the unexpected");
    }
    #[test]
    #[should_panic(expected = "Maximum extensibility reached.")]
    fn calculate_helmholtz_free_energy_density() {
        Gent::new(GENTPARAMETERS)
            .calculate_helmholtz_free_energy_density(&DeformationGradient::new([
                [16.0, 0.00, 0.00],
                [0.0, 0.25, 0.00],
                [0.0, 0.00, 0.25],
            ]))
            .expect("the unexpected");
    }
}
