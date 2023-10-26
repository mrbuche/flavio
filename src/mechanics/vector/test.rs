use super::
{
    Vector,
    VectorTrait
};

#[test]
fn todo()
{
    todo!()
    // &Vector::<8, 0>::zero() * &Vector::<8, 1>::zero();
    // how to test enforcement that contractions over different configurations are forbidden?
}

#[test]
fn zero()
{
    Vector::<8, 1>::zero().iter()
    .for_each(|vector_i|
        assert_eq!(vector_i, &0.0)
    );
}