use super::
{
    AlmansiHamel,
    Elastic,
    super::test::
    {
        ALMANSIHAMELPARAMETERS,
        test_elastic_constitutive_model,
        test_elastic_only_constitutive_model_constructed
    }
};

test_elastic_constitutive_model!(
    AlmansiHamel,
    ALMANSIHAMELPARAMETERS,
    AlmansiHamel::new(ALMANSIHAMELPARAMETERS)
);

test_elastic_only_constitutive_model_constructed!(
    AlmansiHamel::new(ALMANSIHAMELPARAMETERS)
);