use crate::constitutive::
{
    Constitutive,
    hybrid::
    {
        Additive,
        Hybrid
    },
    solid::
    {
        hyperelastic::
        {
            Fung,
            Gent,
            test::
            {
                FUNGPARAMETERS,
                GENTPARAMETERS
            }
        }
    }
};

#[test]
fn todo()
{
    let _ = Additive::construct(
        Fung::new(FUNGPARAMETERS),
        Gent::new(GENTPARAMETERS)
    );
    // send through hyperelastic tests macros for each hyperelastic combo?
    // and hydrid speicific tests, like gets()?
    // and moduli if you can get them working?
    // put hybrid-specific testing macros (or additive, multiplicative) in hybrid/tests
    // and then for hyperelastic in one down
    // then apply the tests another one down
}
