#[cfg(test)]
mod test;

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

pub fn inverse_langevin_approximate(y: f64) -> f64
{
    (2.14234*y.powi(3) - 4.22785*y.powi(2) + 3.0*y)/(1.0 - y)/(0.71716*y.powi(3) - 0.41103*y.powi(2) - 0.39165*y + 1.0)
}

pub fn langevin(x: f64) -> f64
{
    1.0/x.tanh() - 1.0/x
}

pub fn langevin_derivative(x: f64) -> f64
{
    1.0/x.powi(2) - 1.0/x.sinh().powi(2)
}