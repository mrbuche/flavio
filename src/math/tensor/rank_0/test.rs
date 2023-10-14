use super::TensorRank0;

#[test]
fn copy()
{
    let a: TensorRank0 = 1.0;
    let b = a;
    assert_eq!(a, b);
}

#[test]
fn add()
{
    let a: TensorRank0 = 1.0;
    let b: TensorRank0 = 2.0;
    assert_eq!(a + b, 3.0);
}

#[test]
fn subtract()
{
    let a: TensorRank0 = 1.0;
    let b: TensorRank0 = 2.0;
    assert_eq!(a - b, -1.0);
}

#[test]
fn multiply()
{
    let a: TensorRank0 = 1.0;
    let b: TensorRank0 = 2.0;
    assert_eq!(a * b, 2.0);
}

#[test]
fn divide()
{
    let a: TensorRank0 = 1.0;
    let b: TensorRank0 = 2.0;
    assert_eq!(a / b, 0.5);
}
