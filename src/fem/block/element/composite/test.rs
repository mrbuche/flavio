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
            },
            test::assert_eq_within_tols
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
        fn partition_of_unity<'a>()
        {
            $element::<AlmansiHamel>::calculate_standard_gradient_operators().iter()
            .for_each(|standard_gradient_operator|{
                let mut sum = [0.00_f64; 3];
                standard_gradient_operator.iter()
                .for_each(|row|
                    row.iter()
                    .zip(sum.iter_mut())
                    .for_each(|(entry, sum_i)|
                        *sum_i += entry
                    )
                );
                sum.iter()
                .for_each(|sum_i|
                    assert_eq_within_tols(sum_i, &0.00)
                )
            })
        }
        #[test]
        fn todo()
        {
            todo!("can probably test sum of shape functions at integration points is unity as well")
        }
        #[test]
        fn todo1()
        {
            todo!("also test that sum of shape function integrals inverse is that 4x4 matrix")
        }
        #[test]
        fn todo2()
        {
            todo!("test composite jacobian somehow?")
        }
        #[test]
        fn todo3()
        {
            todo!("and what about subtet jacobians, parameteric operators")
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
