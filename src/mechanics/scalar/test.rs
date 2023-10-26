use super::Scalar;

#[test]
fn tensor_rank_0()
{
    let _: Scalar;
}

#[test]
fn copy()
{
    let a: Scalar = 1.0;
    let b = a;
    assert_eq!(a, b);
}

#[test]
fn add()
{
    let a: Scalar = 1.0;
    let b: Scalar = 2.0;
    assert_eq!(a + b, 3.0);
}

#[test]
fn subtract()
{
    let a: Scalar = 1.0;
    let b: Scalar = 2.0;
    assert_eq!(a - b, -1.0);
}

#[test]
fn multiply()
{
    let a: Scalar = 1.0;
    let b: Scalar = 2.0;
    assert_eq!(a * b, 2.0);
}

#[test]
fn divide()
{
    let a: Scalar = 1.0;
    let b: Scalar = 2.0;
    assert_eq!(a / b, 0.5);
}
