use super::
{
    *, super::test::*
};

test_elastic_constitutive_model_nu!(
    AlmansiHamel,
    ALMANSIHAMELPARAMETERS,
    AlmansiHamel::new(ALMANSIHAMELPARAMETERS)
);

test_elastic_constitutive_model!(
    AlmansiHamel,
    ALMANSIHAMELPARAMETERS,
    AlmansiHamel::new(ALMANSIHAMELPARAMETERS)
);

test_elastic_only_constitutive_model_constructed!(
    AlmansiHamel::new(ALMANSIHAMELPARAMETERS)
);