#[cfg(test)]
mod test;

use crate::
{
    ABS_TOL,
    constitutive::
    {
        Constitutive,
        Parameters,
        hybrid::
        {
            Hybrid,
            Multiplicative,
            MultiplicativeTrait
        },
        solid::
        {
            Solid,
            elastic::Elastic
        }
    },
    math::
    {
        TensorRank2,
        TensorRank2Trait
    },
    mechanics::
    {
        CauchyStress,
        CauchyTangentStiffness,
        DeformationGradient,
        FirstPiolaKirchoffStress,
        FirstPiolaKirchoffTangentStiffness,
        Scalar,
        SecondPiolaKirchoffStress,
        SecondPiolaKirchoffTangentStiffness
    }
};

/// Constitutive model implementation of hybrid elastic constitutive models that are based on the multiplicative decomposition.
impl<'a, C1: Elastic<'a>, C2: Elastic<'a>> Constitutive<'a> for Multiplicative<C1, C2>
{
    /// Dummy method that will panic, use [Self::construct()] instead.
    fn new(_parameters: Parameters<'a>) -> Self
    {
        panic!()
    }
}

/// Solid constitutive model implementation of hybrid elastic constitutive models that are based on the multiplicative decomposition.
impl<'a, C1: Elastic<'a>, C2: Elastic<'a>> Solid<'a> for Multiplicative<C1, C2>
{
    /// Dummy method that will panic.
    fn get_bulk_modulus(&self) -> &Scalar
    {
        panic!()
    }
    /// Dummy method that will panic.
    fn get_shear_modulus(&self) -> &Scalar
    {
        panic!()
    }
}

/// Elastic constitutive model implementation of hybrid elastic constitutive models that are based on the multiplicative decomposition.
impl<'a, C1: Elastic<'a>, C2: Elastic<'a>> Elastic<'a> for Multiplicative<C1, C2>
{
    /// Calculates and returns the Cauchy stress.
    ///
    /// ```math
    /// \boldsymbol{\sigma}(\mathbf{F}) = \frac{1}{J_2}\,\boldsymbol{\sigma}_1(\mathbf{F}_1)
    /// ```
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        let (deformation_gradient_1, deformation_gradient_2) = self.calculate_deformation_gradients(deformation_gradient);
        self.get_constitutive_model_1().calculate_cauchy_stress(&deformation_gradient_1) / deformation_gradient_2.determinant()
    }
    /// Dummy method that will panic.
    fn calculate_cauchy_tangent_stiffness(&self, _: &DeformationGradient) -> CauchyTangentStiffness
    {
        panic!()
    }
    /// Calculates and returns the first Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathbf{P}(\mathbf{F}) = \mathbf{P}_1(\mathbf{F}_1)\cdot\mathbf{F}_2^{-T}
    /// ```
    fn calculate_first_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient) -> FirstPiolaKirchoffStress
    {
        let (deformation_gradient_1, deformation_gradient_2) = self.calculate_deformation_gradients(deformation_gradient);
        let deformation_gradient_2_inverse_transpose: TensorRank2<3, 0, 0> = deformation_gradient_2.inverse_transpose().into();
        self.get_constitutive_model_1().calculate_first_piola_kirchoff_stress(&deformation_gradient_1) * deformation_gradient_2_inverse_transpose
    }
    /// Dummy method that will panic.
    fn calculate_first_piola_kirchoff_tangent_stiffness(&self, _: &DeformationGradient) -> FirstPiolaKirchoffTangentStiffness
    {
        panic!()
    }
    /// Calculates and returns the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathbf{S}(\mathbf{F}) = \mathbf{F}_2^{-1}\cdot\mathbf{S}_1(\mathbf{F}_1)\cdot\mathbf{F}_2^{-T}
    /// ```
    fn calculate_second_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient) -> SecondPiolaKirchoffStress
    {
        let (deformation_gradient_1, deformation_gradient_2) = self.calculate_deformation_gradients(deformation_gradient);
        let deformation_gradient_2_inverse: TensorRank2<3, 0, 0> = deformation_gradient_2.inverse().into();
        &deformation_gradient_2_inverse * self.get_constitutive_model_1().calculate_second_piola_kirchoff_stress(&deformation_gradient_1) * deformation_gradient_2_inverse.transpose()
    }
    /// Dummy method that will panic.
    fn calculate_second_piola_kirchoff_tangent_stiffness(&self, _: &DeformationGradient) -> SecondPiolaKirchoffTangentStiffness
    {
        panic!()
    }
}

/// Multiplicative hybrid constitutive model implementation of hybrid elastic constitutive models that are based on the multiplicative decomposition.
impl<'a, C1: Elastic<'a>, C2: Elastic<'a>> MultiplicativeTrait for Multiplicative<C1, C2>
{
    fn calculate_deformation_gradients(&self, deformation_gradient: &DeformationGradient) -> (DeformationGradient, DeformationGradient)
    {
        if deformation_gradient.is_identity()
        {
            (DeformationGradient::identity(), DeformationGradient::identity())
        }
        else
        {
            let mut deformation_gradient_1 = DeformationGradient::identity();
            let mut deformation_gradient_2 = deformation_gradient * 1.0;
            let mut deformation_gradient_2_old = DeformationGradient::identity();
            let mut deformation_gradient_2_inverse_transpose: TensorRank2<3, 0, 0>;
            let mut residual: FirstPiolaKirchoffStress;
            let mut residual_increment: FirstPiolaKirchoffStress;
            let mut residual_norm = 1.0;
            let mut residual_old = FirstPiolaKirchoffStress::zero();
            let mut right_hand_side: FirstPiolaKirchoffStress;
            let mut steps: u8 = 0;
            let steps_maximum: u8 = 50;
            let mut step_size: Scalar;
            while steps < steps_maximum
            {
                deformation_gradient_1 = (deformation_gradient * deformation_gradient_2.inverse()).into();
                deformation_gradient_2_inverse_transpose = deformation_gradient_2.inverse_transpose().into();
                right_hand_side = (deformation_gradient_1.transpose() * self.get_constitutive_model_1().calculate_first_piola_kirchoff_stress(&deformation_gradient_1) * deformation_gradient_2_inverse_transpose).into();
                residual = self.get_constitutive_model_2().calculate_first_piola_kirchoff_stress(&deformation_gradient_2) - right_hand_side;
                residual_norm = residual.norm();
                if residual_norm >= ABS_TOL
                {
                    residual_increment = residual_old - &residual;
                    step_size = (deformation_gradient_2_old - &deformation_gradient_2).full_contraction(&residual_increment).abs() / residual_increment.norm_squared();
                    deformation_gradient_2_old = &deformation_gradient_2 * 1.0;
                    residual_old = &residual * 1.0;
                    deformation_gradient_2 -= residual * step_size;
                }
                else
                {
                    break
                }
                steps += 1;
            }
            if residual_norm >= ABS_TOL && steps == steps_maximum
            {
                panic!("The maximum number of steps was reached before the tolerance was satisfied.")
            }
            (deformation_gradient_1, deformation_gradient_2)
        }
    }
}
