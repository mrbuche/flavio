use crate::test::assert_eq_within_tols;
use super::*;

mod inverse_langevin
{
    use super::*;
    #[test]
    fn range()
    {
        let mut gamma = -1.0;
        let length = 9999;
        let d_gamma = 2.0/((length + 1) as f64);
        for _ in 0..length
        {
            gamma += d_gamma;
            assert_eq_within_tols(&langevin(inverse_langevin(gamma)), &gamma)
        }
    }
    #[test]
    fn zero()
    {
        assert_eq!(inverse_langevin(0.0), 0.0)
    }
}
