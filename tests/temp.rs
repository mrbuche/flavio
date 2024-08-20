#[cfg(feature = "constitutive")]
mod temporary {
    use flavio::constitutive::{
        solid::{elastic::Elastic, hyperelastic::NeoHookean},
        Constitutive,
    };
    #[test]
    fn neo_hookean() {
        let model = NeoHookean::new(&[13.0, 3.0]);
        let stress = model.solve_uniaxial_tension(&13.0).expect("the unexpected");
        println!("{}", stress)
    }
}
