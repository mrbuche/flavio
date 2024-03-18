use crate::constitutive::solid::hyperelastic::
{
    NeoHookean,
    test::NEOHOOKEANPARAMETERS
};
use crate::test::assert_eq_within_tols;
use super::*;

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.0, 0.0, 0.5],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5]
    ])
}

fn get_element<'a>() -> Tetrahedron::<NeoHookean<'a>>
{
    Tetrahedron::new(
        NEOHOOKEANPARAMETERS, get_reference_coordinates()
    )
}

#[test]
fn partition_of_unity<'a>()
{
    Tetrahedron::<NeoHookean<'a>>::calculate_standard_gradient_operators()
    .iter()
    .for_each(|standard_gradient_operator|{
        let mut sum = [0.0_f64; 3];
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
            assert_eq_within_tols(sum_i, &0.0)
        )
    })
}

#[test]
fn todo()
{
    todo!()
}