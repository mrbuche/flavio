#[cfg(feature = "mechanics")]
mod temporary
{
    use flavio::
    {
        constitutive::solid::hyperelastic::NeoHookean,
        fem::*,
        math::*,
        mechanics::*
    };
    type NodalCoordinates<const D: usize> = CurrentCoordinates<D>;
    type ReferenceNodalCoordinates<const D: usize> = ReferenceCoordinates<D>;
    fn get_reference_coordinates() -> ReferenceNodalCoordinates<12>
    {
        ReferenceNodalCoordinates::new([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.75_f64.sqrt(), 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.75_f64.sqrt(), 0.0],
            [0.5, 0.0, 0.0],
            [0.75, 0.5 * 0.75_f64.sqrt(), 0.0],
            [0.25, 0.5 * 0.75_f64.sqrt(), 0.0],
            [0.5, 0.0, 0.0],
            [0.75, 0.5 * 0.75_f64.sqrt(), 0.0],
            [0.25, 0.5 * 0.75_f64.sqrt(), 0.0]
        ])
    }
    fn get_deformed_coordinates() -> NodalCoordinates<12>
    {
        NodalCoordinates::new([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.75_f64.sqrt(), 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.5, 0.75_f64.sqrt(), 1.0],
            [0.5, 0.0, 0.0],
            [0.75, 0.5 * 0.75_f64.sqrt(), 0.0],
            [0.25, 0.5 * 0.75_f64.sqrt(), 0.0],
            [0.5, 0.0, 1.0],
            [0.75, 0.5 * 0.75_f64.sqrt(), 1.0],
            [0.25, 0.5 * 0.75_f64.sqrt(), 1.0]
        ])
    }
    #[test]
    fn temporary()
    {
        let element = CompositeWedgeLocalization::<NeoHookean>::new(
            &[13.0, 3.0],
            get_reference_coordinates(),
            &1.23
        );
        element.calculate_nodal_forces(
            &get_deformed_coordinates()
        ).iter().for_each(|nodal_force|
            println!("{:?}", (nodal_force[0], nodal_force[1], nodal_force[2]))
        );
    }
}
