macro_rules! test_linear_localization_element
{
    ($element: ident) =>
    {
        crate::fem::block::element::test::setup_for_localization_elements!($element);
        crate::fem::block::element::linear::test::test_linear_element_inner!($element);
        crate::fem::block::element::linear::surface::test::test_linear_surface_element_inner!($element);
    }
}
pub(crate) use test_linear_localization_element;

#[test]
fn test_localization_specifics()
{
    todo!()
}