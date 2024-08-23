use super::super::test::*;
use super::*;

use_elastic_macros!();

test_solid_hyperelastic_constitutive_model!(
    ArrudaBoyce,
    ARRUDABOYCEPARAMETERS,
    ArrudaBoyce::new(ARRUDABOYCEPARAMETERS)
);

test_solve_uniaxial_tension!(ArrudaBoyce::new(ARRUDABOYCEPARAMETERS));

#[test]
fn get_number_of_links() {
    assert_eq!(
        &ARRUDABOYCEPARAMETERS[2],
        ArrudaBoyce::new(ARRUDABOYCEPARAMETERS).get_number_of_links()
    )
}

mod maximum_extensibility {
    use super::*;
    #[test]
    #[should_panic(expected = "Maximum extensibility reached.")]
    fn calculate_cauchy_stress() {
        ArrudaBoyce::new(ARRUDABOYCEPARAMETERS)
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
        ArrudaBoyce::new(ARRUDABOYCEPARAMETERS)
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
        ArrudaBoyce::new(ARRUDABOYCEPARAMETERS)
            .calculate_helmholtz_free_energy_density(&DeformationGradient::new([
                [16.0, 0.00, 0.00],
                [0.0, 0.25, 0.00],
                [0.0, 0.00, 0.25],
            ]))
            .expect("the unexpected");
    }
}
