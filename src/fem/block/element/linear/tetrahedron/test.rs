use crate::fem::block::element::
{
    test::test_finite_element,
    linear::
    {
        test::test_linear_finite_element
    }
};
use super::*;

fn get_element_standard_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0]
    ])
}

test_finite_element!(LinearTetrahedron);
test_linear_finite_element!(LinearTetrahedron);

#[test]
fn specific_tests_for_linear_tets()
{
    todo!()
}