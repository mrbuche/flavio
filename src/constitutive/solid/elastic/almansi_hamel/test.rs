use super::
{
    *, super::test::*
};

test_solid_constitutive_model!(
    AlmansiHamel,
    ALMANSIHAMELPARAMETERS,
    AlmansiHamel::new(ALMANSIHAMELPARAMETERS)
);

test_solid_elastic_constitutive_model!(
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