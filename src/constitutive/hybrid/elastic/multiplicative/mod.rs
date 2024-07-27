#[cfg(test)]
mod test;

use crate::
{
    constitutive::
    {
        Constitutive,
        Parameters,
        hybrid::
        {
            Multiplicative,
            MultiplicativeTrait
        },
        solid::
        {
            Solid,
            elastic::Elastic
        }
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

// Inherent implementation of hybrid elastic constitutive models that are based on the multiplicative decomposition.
// impl<'a, C1: Elastic<'a>, C2: Elastic<'a>> Multiplicative<C1, C2>
// {
//     /// Calculates and returns the first deformation gradient.
//     ///
//     /// ```math
//     /// ???
//     /// ```
//     fn calculate_deformation_gradient_1(deformation_gradient: &DeformationGradient) -> DeformationGradient
//     {
//         todo!("Do minimization for F1 or F2 based on what is better for use below; update description and equation.
//                Note that this should actually be an Elastic trait method,
//                and below should be a Hyperelastic trait method.")
//     }
//     fn calculate_deformation_gradients(deformation_gradient: &DeformationGradient) -> (DeformationGradient, DeformationGradient)
//     {
//         todo!("Need both for Helmholtz free energy; avoid code duplication with above.")
//     }
// }

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
    fn get_bulk_modulus(&self) -> &Scalar
    {
        panic!()
    }
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
    /// \boldsymbol{\sigma}(\mathbf{F}) = ?
    /// ```
    fn calculate_cauchy_stress(&self, deformation_gradient: &DeformationGradient) -> CauchyStress
    {
        todo!()
    }
    /// Calculates and returns the tangent stiffness associated with the Cauchy stress.
    ///
    /// ```math
    /// \mathcal{T}(\mathbf{F}) = ?
    /// ```
    fn calculate_cauchy_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> CauchyTangentStiffness
    {
        todo!()
    }
    /// Calculates and returns the first Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathbf{P}(\mathbf{F}) = ?
    /// ```
    fn calculate_first_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient) -> FirstPiolaKirchoffStress
    {
        todo!()
    }
    /// Calculates and returns the tangent stiffness associated with the first Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathcal{C}(\mathbf{F}) = ?
    /// ```
    fn calculate_first_piola_kirchoff_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> FirstPiolaKirchoffTangentStiffness
    {
        todo!()
    }
    /// Calculates and returns the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathbf{S}(\mathbf{F}) = ?
    /// ```
    fn calculate_second_piola_kirchoff_stress(&self, deformation_gradient: &DeformationGradient) -> SecondPiolaKirchoffStress
    {
        todo!()
    }
    /// Calculates and returns the tangent stiffness associated with the second Piola-Kirchoff stress.
    ///
    /// ```math
    /// \mathcal{G}(\mathbf{F}) = ?
    /// ```
    fn calculate_second_piola_kirchoff_tangent_stiffness(&self, deformation_gradient: &DeformationGradient) -> SecondPiolaKirchoffTangentStiffness
    {
        todo!()
    }
}

/// Multiplicative hybrid constitutive model implementation of hybrid elastic constitutive models that are based on the multiplicative decomposition.
impl<'a, C1: Elastic<'a>, C2: Elastic<'a>> MultiplicativeTrait for Multiplicative<C1, C2>
{}
