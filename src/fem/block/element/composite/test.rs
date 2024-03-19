macro_rules! test_composite_element
{
    ($element: ident) =>
    {
        crate::fem::block::element::composite::test::setup_for_composite_elements!($element);
        crate::fem::block::element::test::test_finite_element!($element);
    }
}
pub(crate) use test_composite_element;

macro_rules! setup_for_composite_elements
{
    ($element: ident) =>
    {
        use crate::
        {
            constitutive::solid::elastic::AlmansiHamel,
            mechanics::test::
            {
                get_deformation_gradient,
                get_deformation_gradient_rate
            }
        };
        fn get_coordinates() -> NodalCoordinates<N>
        {
            get_deformation_gradient() * get_reference_coordinates()
        }
        fn get_velocities() -> NodalVelocities<N>
        {
            get_deformation_gradient_rate() * get_reference_coordinates()
        }
        #[test]
        fn size()
        {
            assert_eq!(
                std::mem::size_of::<$element::<AlmansiHamel>>(),
                std::mem::size_of::<[AlmansiHamel; G]>()
                + std::mem::size_of::<ProjectedGradientVectors<G, N>>()
            )
        }
        crate::fem::block::element::test::setup_for_element_tests_any_element!($element);
    }
}
pub(crate) use setup_for_composite_elements;

#[test]
fn todo()
{
    todo!()
}