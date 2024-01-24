#![cfg(feature = "doc")]

use flavio::
{
    constitutive::
    {
        ConstitutiveModel,
        hyperelastic::
        {
            ArrudaBoyceModel,
            FungModel,
            GentModel,
            MooneyRivlinModel,
            NeoHookeanModel,
            SaintVenantKirchoffModel,
            YeohModel
        }
    },
    math::TensorRank2Trait,
    mechanics::DeformationGradient
};

use plotpy::
{
    linspace,
    Curve,
    Plot
};

fn incompressible_uniaxial_stretch<'a, C: ConstitutiveModel<'a>>(lambda: &Vec<f64>, model: C) -> Vec<f64>
{
    lambda.iter().map(|lambda_i|
        model.calculate_cauchy_stress(
            &DeformationGradient::new([
                [*lambda_i, 0.0, 0.0],
                [0.0, 1.0/lambda_i.sqrt(), 0.0],
                [0.0, 0.0, 1.0/lambda_i.sqrt()],
            ])
        )[0][0]
    ).collect()
}

#[test]
fn plots()
{
    let mut plot = Plot::new();
    let lambda = linspace(0.1, 5.0, 333);
    let mut curve1 = Curve::new();
    curve1.set_label("Arruda-Boyce").set_line_width(2.0);
    curve1.draw(&lambda, &incompressible_uniaxial_stretch(&lambda, ArrudaBoyceModel::new(&[1.0, 1.0, 16.0])));
    let mut curve2 = Curve::new();
    curve2.set_label("Fung").set_line_width(2.0);
    curve2.draw(&lambda, &incompressible_uniaxial_stretch(&lambda, FungModel::new(&[1.0, 1.0, 1.0, 0.1])));
    let mut curve3 = Curve::new();
    curve3.set_label("Gent").set_line_width(2.0);
    curve3.draw(&lambda, &incompressible_uniaxial_stretch(&lambda, GentModel::new(&[1.0, 1.0, 30.0])));
    let mut curve4 = Curve::new();
    curve4.set_label("Mooney-Rivlin").set_line_width(2.0);
    curve4.draw(&lambda, &incompressible_uniaxial_stretch(&lambda, MooneyRivlinModel::new(&[1.0, 1.0, 0.3])));
    let mut curve5 = Curve::new();
    curve5.set_label("Neo-Hookean").set_line_width(2.0);
    curve5.draw(&lambda, &incompressible_uniaxial_stretch(&lambda, NeoHookeanModel::new(&[1.0, 1.0])));
    let mut curve6 = Curve::new();
    curve6.set_label("St. Venant-Kirchoff").set_line_width(2.0);
    curve6.draw(&lambda, &incompressible_uniaxial_stretch(&lambda, SaintVenantKirchoffModel::new(&[1.0, 1.0])));
    let mut curve7 = Curve::new();
    curve7.set_label("Yeoh").set_line_width(2.0);
    curve7.draw(&lambda, &incompressible_uniaxial_stretch(&lambda, YeohModel::new(&[1.0, 1.0, 0.1, 0.02])));
    plot.add(&curve1).add(&curve2).add(&curve3).add(&curve4).add(&curve5).add(&curve6).add(&curve7);
    plot.grid_labels_legend("λ", "σ/μ").set_xmin(0.0).set_xmax(5.0).set_ymin(-10.0).set_ymax(20.0);
    let _ = plot.save("target/doc/flavio/constitutive/hyperelastic/models.svg");
}