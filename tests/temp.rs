#[cfg(feature = "fem")]
mod temporary {
    use flavio::{
        constitutive::solid::hyperelastic::NeoHookean,
        fem::{FiniteElement, HyperelasticFiniteElement, LinearTetrahedron},
        math::{TensorRank1ListTrait, TensorRank2Trait},
        mechanics::{CurrentCoordinates, DeformationGradient, ReferenceCoordinates},
    };
    fn get_coordinates() -> CurrentCoordinates<4> {
        let mut deformation_gradient = DeformationGradient::identity();
        deformation_gradient[0][0] *= -1.0;
        deformation_gradient * get_reference_coordinates()
    }
    fn get_reference_coordinates() -> ReferenceCoordinates<4> {
        ReferenceCoordinates::new([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ])
    }
    #[test]
    fn invalid_jacobian() {
        let element =
            LinearTetrahedron::<NeoHookean>::new(&[13.0, 3.0], get_reference_coordinates());
        element
            .calculate_helmholtz_free_energy(&get_coordinates())
            .expect("the unexpected");
    }
}
