#[test]
fn size()
{
    assert_eq!(
        std::mem::size_of::<super::ConstitutiveModelParameters>(), 16
    )
}