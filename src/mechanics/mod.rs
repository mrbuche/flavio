#[cfg(test)]
pub mod test;

use crate::math::
{
    TensorRank0,
    TensorRank1,
    TensorRank1List,
    TensorRank2,
    TensorRank2List,
    TensorRank4
};

/// The Cauchy stress $`\boldsymbol{\sigma}`$.
pub type CauchyStress = TensorRank2<3, 1, 1>;

/// The tangent stiffness associated with the Cauchy stress $`\boldsymbol{\mathcal{T}}`$.
pub type CauchyTangentStiffness = TensorRank4<3, 1, 1, 1, 0>;

/// A coordinate in the current configuration.
pub type CurrentCoordinate = TensorRank1<3, 0>;

/// A list of coordinates in the current configuration.
pub type CurrentCoordinates<const L: usize> = TensorRank1List<3, 0, L>;

/// The deformation gradient $`\mathbf{F}`$.
pub type DeformationGradient = TensorRank2<3, 1, 0>;

/// A list of deformation gradients.
pub type DeformationGradients<const L: usize> = TensorRank2List<3, 1, 0, L>;

/// The first Piola-Kirchoff stress $`\mathbf{P}`$.
pub type FirstPiolaKirchoffStress = TensorRank2<3, 1, 0>;

/// A list of first Piola-Kirchoff stresses.
pub type FirstPiolaKirchoffStresses<const L: usize> = TensorRank2List<3, 1, 0, L>;

/// The tangent stiffness associated with the first Piola-Kirchoff stress $`\boldsymbol{\mathcal{C}}`$.
pub type FirstPiolaKirchoffTangentStiffness = TensorRank4<3, 1, 0, 1, 0>;

/// A force.
pub type Force = TensorRank1<3, 1>;

/// A list of forces.
pub type Forces<const L: usize> = TensorRank1List<3, 1, L>;

/// The left Cauchy-Green deformation $`\mathbf{B}`$.
pub type LeftCauchyGreenDeformation = TensorRank2<3, 1, 1>;

/// A coordinate in the reference configuration.
pub type ReferenceCoordinate = TensorRank1<3, 0>;

/// A list of coordinates in the reference configuration.
pub type ReferenceCoordinates<const L: usize> = TensorRank1List<3, 0, L>;

/// The right Cauchy-Green deformation $`\mathbf{C}`$.
pub type RightCauchyGreenDeformation = TensorRank2<3, 0, 0>;

/// The rotation of the current configuration $`\mathbf{Q}`$.
pub type RotationCurrentConfiguration = TensorRank2<3, 1, 1>;

/// The rotation of the reference configuration $`\mathbf{Q}_0`$.
pub type RotationReferenceConfiguration = TensorRank2<3, 0, 0>;

/// An arbitrary scalar.
pub type Scalar = TensorRank0;

/// A stiffness resulting from a force.
pub type Stiffness = TensorRank2<3, 1, 1>;

/// A list of stiffnesses.
pub type Stiffnesses<const L: usize> = TensorRank2List<3, 1, 1, L>;
