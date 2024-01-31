use super::
{
    AlmansiHamelModel,
    ElasticConstitutiveModel,
    super::test::
    {
        ALMANSIHAMELPARAMETERS,
        test_elastic_constitutive_model,
        test_elastic_only_constitutive_model_constructed
    }
};

test_elastic_constitutive_model!(
    AlmansiHamelModel,
    ALMANSIHAMELPARAMETERS,
    AlmansiHamelModel::new(ALMANSIHAMELPARAMETERS)
);

test_elastic_only_constitutive_model_constructed!(
    AlmansiHamelModel::new(ALMANSIHAMELPARAMETERS)
);