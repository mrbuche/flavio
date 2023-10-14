use super::TensorRank1;

#[test]
fn copy()
{
    let a = TensorRank1([1.0, 2.0, 3.0]);
    let b = a;
    a.iter().zip(b.iter()).for_each(|(a_i, b_i)|
        assert_eq!(a_i, b_i)
    );
}

#[test]
fn todo()
{
    todo!("Make comprehensive tests, and traits for <D> and <3> instead of base impls in order to use for mechanics.");
}
