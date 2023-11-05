#[cfg(test)]
mod test;

/// Returns the inverse Langevin function.
///
/// Improves [inverse_langevin_approximate()] using two iterations of Newton's method.
pub fn inverse_langevin(y: f64) -> f64
{
    let y_abs = y.abs();
    if y_abs >= 1.0
    {
        panic!()
    }
    else if y_abs <= 1e-3
    {
        3.0 * y + 9.0/5.0 * y.powi(3)
    }
    else
    {
        let mut x = inverse_langevin_approximate(y_abs);
        for _ in 0..2
        {
            x += (y_abs - langevin(x)) / langevin_derivative(x);
        }
        if y < 0.0
        {
            -x
        }
        else
        {
            x
        }
    }
}

/// Returns an approximation of the inverse Langevin function.
///
/// The approximation given by [Jedynak](https://doi.org/10.1177%2F1081286518811395) has a maximum relative error of 0.082%.
///
/// ```math
/// \mathcal{L}^{-1}(y) \approx \frac{2.14234 y^3 - 4.22785 y^2 + 3y}{(1 - y)(0.71716 y^3 - 0.41103 y^2 - 0.39165 y + 1)}
/// ```
pub fn inverse_langevin_approximate(y: f64) -> f64
{
    (2.14234*y.powi(3) - 4.22785*y.powi(2) + 3.0*y)/(1.0 - y)/(0.71716*y.powi(3) - 0.41103*y.powi(2) - 0.39165*y + 1.0)
}

/// Returns the Langevin function.
///
/// ```math
/// \mathcal{L}(x) = \coth(x) - x^{-1}
/// ```
pub fn langevin(x: f64) -> f64
{
    1.0/x.tanh() - 1.0/x
}

/// Returns derivative of the Langevin function.
///
/// ```math
/// \mathcal{L}'(x) = x^{-2} - \sinh^{-2}(x)
/// ```
pub fn langevin_derivative(x: f64) -> f64
{
    1.0/x.powi(2) - 1.0/x.sinh().powi(2)
}