use crate::fem::block::element::linear::surface::test::test_linear_surface_element;
use super::*;

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]
    ])
}

test_linear_surface_element!(Triangle);

fn get_normal_rate_from_finite_difference() -> Normal<1>
{
    let mut finite_difference = 0.0;
    (0..3).map(|i|
        get_velocities().iter().enumerate()
        .map(|(a, velocity_a)|
            velocity_a.iter().enumerate()
            .map(|(k, velocity_a_k)|{
                let mut coordinates = get_coordinates();
                coordinates[a][k] += 0.5 * EPSILON;
                finite_difference = Triangle::<AlmansiHamel>::calculate_normal(
                    &coordinates
                )[i];
                coordinates[a][k] -= EPSILON;
                finite_difference -= Triangle::<AlmansiHamel>::calculate_normal(
                    &coordinates
                )[i];
                finite_difference/EPSILON * velocity_a_k
            }).sum::<Scalar>()
        ).sum()
    ).collect()
}

fn get_normal_rate(is_deformed: bool) -> NormalRate
{
    if is_deformed
    {
        Triangle::<AlmansiHamel>::calculate_normal_rate(
            &get_coordinates(),
            &get_velocities()
        )
    }
    else
    {
        Triangle::<AlmansiHamel>::calculate_normal_rate(
            &get_reference_coordinates().convert(),
            &NodalVelocities::zero()
        )
    }
}