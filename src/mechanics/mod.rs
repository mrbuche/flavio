#[cfg(test)]
pub mod test;

use crate::math::
{
    TensorRank0,
    TensorRank2,
    TensorRank4
};

/// The Cauchy stress $`\boldsymbol{\sigma}`$.
pub type CauchyStress = TensorRank2<3, 1, 1>;

/// The tangent stiffness associated with the Cauchy stress $`\boldsymbol{\mathcal{T}}`$.
pub type CauchyTangentStiffness = TensorRank4<3, 1, 1, 1, 0>;

/// The deformation gradient $`\mathbf{F}`$.
pub type DeformationGradient = TensorRank2<3, 1, 0>;

/// The first Piola-Kirchoff stress $`\mathbf{P}`$.
pub type FirstPiolaKirchoffStress = TensorRank2<3, 1, 0>;

/// The tangent stiffness associated with the first Piola-Kirchoff stress $`\boldsymbol{\mathcal{C}}`$.
pub type FirstPiolaKirchoffTangentStiffness = TensorRank4<3, 1, 0, 1, 0>;

/// The left Cauchy-Green deformation $`\mathbf{B}`$.
pub type LeftCauchyGreenDeformation = TensorRank2<3, 1, 1>;

/// The right Cauchy-Green deformation $`\mathbf{C}`$.
pub type RightCauchyGreenDeformation = TensorRank2<3, 0, 0>;

/// The rotation of the current configuration $`\mathbf{Q}`$.
pub type RotationCurrentConfiguration = TensorRank2<3, 1, 1>;

/// The rotation of the reference configuration $`\mathbf{Q}_0`$.
pub type RotationReferenceConfiguration = TensorRank2<3, 0, 0>;

/// An arbitrary scalar.
pub type Scalar = TensorRank0;
