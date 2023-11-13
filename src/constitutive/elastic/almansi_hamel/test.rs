use crate::constitutive::
{
    test::ALMANSIHAMELPARAMETERS,
    elastic::test::
    {
        test_elastic_constitutive_model,
        test_elastic_only_constitutive_model_constructed
    }
};
use super::
{
    AlmansiHamelModel,
    ElasticConstitutiveModel
};

test_elastic_constitutive_model!(
    AlmansiHamelModel,
    ALMANSIHAMELPARAMETERS,
    AlmansiHamelModel::new(ALMANSIHAMELPARAMETERS)
);

test_elastic_only_constitutive_model_constructed!(
    AlmansiHamelModel::new(ALMANSIHAMELPARAMETERS)
);