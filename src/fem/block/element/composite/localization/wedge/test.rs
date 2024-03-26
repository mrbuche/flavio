use super::*;

fn get_reference_coordinates() -> ReferenceNodalCoordinates<N>
{
    ReferenceNodalCoordinates::new([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.5, 0.0, 0.0]
    ])
}

#[test]
fn test_composite_localization_element()
{
    todo!()
}

#[test]
fn test_finite_element_block()
{
    todo!()
}

fn get_element<'a>() -> Wedge<crate::constitutive::solid::hyperelastic::NeoHookean<'a>>
{
    Wedge::new(
        crate::constitutive::solid::hyperelastic::test::NEOHOOKEANPARAMETERS,
        get_reference_coordinates()
    )
}

#[test]
fn temporary_1()
{
    get_element().calculate_deformation_gradients(
        &get_reference_coordinates().convert()
    ).iter()
    .for_each(|deformation_gradient|
        deformation_gradient.iter()
        .enumerate()
        .for_each(|(i, deformation_gradient_i)|
            deformation_gradient_i.iter()
            .enumerate()
            .for_each(|(j, deformation_gradient_ij)|
                crate::test::assert_eq_within_tols(
                    deformation_gradient_ij, &(((i == j) as u8) as Scalar)
                )
            )
        )
    )
}

fn get_deformation_gradient_planar() -> crate::mechanics::DeformationGradient
{
    crate::mechanics::DeformationGradient::new([
        [0.9, 0.5, 0.0],
        [0.7, 0.9, 0.0],
        [0.0, 0.0, 1.0]
    ])
}

fn get_coordinates_from_deformation_gradient_planar() -> NodalCoordinates<N>
{
    NodalCoordinates::new([
        [0.90, 0.70, 0.0],
        [0.50, 0.90, 0.0],
        [0.00, 0.00, 0.0],
        [0.90, 0.70, 0.0],
        [0.50, 0.90, 0.0],
        [0.00, 0.00, 0.0],
        [0.70, 0.80, 0.0],
        [0.25, 0.45, 0.0],
        [0.45, 0.35, 0.0],
        [0.70, 0.80, 0.0],
        [0.25, 0.45, 0.0],
        [0.45, 0.35, 0.0]
    ])
}

#[test]
fn temporary_2()
{
    get_element().calculate_deformation_gradients(
        &get_coordinates_from_deformation_gradient_planar()
    ).iter()
    .for_each(|calculated_deformation_gradient|
        calculated_deformation_gradient.iter()
        .zip(get_deformation_gradient_planar().iter())
        .for_each(|(calculated_deformation_gradient_i, deformation_gradient_i)|
            calculated_deformation_gradient_i.iter()
            .zip(deformation_gradient_i.iter())
            .for_each(|(calculated_deformation_gradient_ij, deformation_gradient_ij)|
                crate::test::assert_eq_within_tols(
                    calculated_deformation_gradient_ij, deformation_gradient_ij
                )
            )
        )
    )
}

fn get_deformation_gradient_both() -> crate::mechanics::DeformationGradient
{
    crate::mechanics::DeformationGradient::new([
        [0.9, 0.5, 0.2],
        [0.7, 0.9, 0.4],
        [0.0, 0.0, 1.6]
    ])
}

fn get_coordinates_from_deformation_gradient_both() -> NodalCoordinates<N>
{
    NodalCoordinates::new([
        [ 0.80,  0.50, -0.3],
        [ 0.40,  0.70, -0.3],
        [-0.10, -0.20, -0.3],
        [ 1.00,  0.90,  0.3],
        [ 0.60,  1.10,  0.3],
        [ 0.10,  0.20,  0.3],
        [ 0.60,  0.60, -0.3],
        [ 0.15,  0.25, -0.3],
        [ 0.35,  0.15, -0.3],
        [ 0.80,  1.00,  0.3],
        [ 0.35,  0.65,  0.3],
        [ 0.55,  0.55,  0.3]
    ])
}

#[test]
fn temporary_3()
{
    get_element().calculate_deformation_gradients(
        &get_coordinates_from_deformation_gradient_both()
    ).iter()
    .for_each(|calculated_deformation_gradient|
        calculated_deformation_gradient.iter()
        .zip(get_deformation_gradient_both().iter())
        .for_each(|(calculated_deformation_gradient_i, deformation_gradient_i)|
            calculated_deformation_gradient_i.iter()
            .zip(deformation_gradient_i.iter())
            .for_each(|(calculated_deformation_gradient_ij, deformation_gradient_ij)|
                crate::test::assert_eq_within_tols(
                    calculated_deformation_gradient_ij, deformation_gradient_ij
                )
            )
        )
    )
}