macro_rules! test_linear_surface_element
{
    ($element: ident) =>
    {
        crate::fem::block::element::test::setup_for_surface_elements!($element);
        crate::fem::block::element::linear::test::test_linear_element_inner!($element);
    }
}
pub(crate) use test_linear_surface_element;