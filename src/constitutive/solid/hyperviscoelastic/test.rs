use crate::
{
    constitutive::solid::viscoelastic::test::ALMANSIHAMELPARAMETERS,
    mechanics::Scalar
};
pub const SAINTVENANTKIRCHOFFPARAMETERS: &[Scalar; 4] = &[ALMANSIHAMELPARAMETERS[0], ALMANSIHAMELPARAMETERS[1], ALMANSIHAMELPARAMETERS[2], ALMANSIHAMELPARAMETERS[3]];

macro_rules! test_hyperviscoelastic_constitutive_model
{
    ($hyperviscoelastic_constitutive_model: ident, $hyperviscoelastic_constitutive_model_parameters: expr, $hyperviscoelastic_constitutive_model_constructed: expr) =>
    {
        use crate::constitutive::solid::
        {
            // viscoelastic::
            // {
            //     Viscoelastic,
            //     test::test_elastic_constitutive_model,
            // },
            hyperviscoelastic::
            {
                // Hyperviscoelastic,
                test::test_hyperviscoelastic_constitutive_model_constructed
            }
        };
        // test_elastic_constitutive_model!($hyperelastic_constitutive_model, $hyperelastic_constitutive_model_parameters, $hyperelastic_constitutive_model_constructed);
        test_hyperviscoelastic_constitutive_model_constructed!($hyperviscoelastic_constitutive_model_constructed);
    }
}
pub(crate) use test_hyperviscoelastic_constitutive_model;

macro_rules! test_hyperviscoelastic_constitutive_model_constructed
{
    ($hyperviscoelastic_constitutive_model_constructed: expr) =>
    {
        mod hyperviscoelastic_only
        {
            mod helmholtz_free_energy_density
            {
                mod deformed
                {
                    #[test]
                    fn finite_difference()
                    {
                        todo!()
                    }
                    #[test]
                    fn minimized()
                    {
                        todo!()
                    }
                    #[test]
                    fn objectivity()
                    {
                        todo!()
                    }
                    #[test]
                    fn positive()
                    {
                        todo!()
                    }
                }
                mod undeformed
                {
                    #[test]
                    fn finite_difference()
                    {
                        todo!()
                    }
                    #[test]
                    fn minimized()
                    {
                        todo!()
                    }
                    #[test]
                    fn zero()
                    {
                        todo!()
                    }
                }
            }
            mod viscous_dissipation
            {
                #[test]
                fn todo()
                {
                    todo!()
                }
            }
        }
    }
}
pub(crate) use test_hyperviscoelastic_constitutive_model_constructed;