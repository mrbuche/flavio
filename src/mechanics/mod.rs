//! Mechanics library.

#[cfg(test)]
pub mod test;

use crate::math::
{
    TensorRank0,
    TensorRank0List,
    TensorRank1,
    TensorRank1List,
    TensorRank1List2D,
    TensorRank2,
    TensorRank2List,
    TensorRank2List2D,
    TensorRank4,
    TensorRank4List
};

/// The Cauchy stress $`\boldsymbol{\sigma}`$.
pub type CauchyStress = TensorRank2<3, 1, 1>;

/// The tangent stiffness associated with the Cauchy stress $`\boldsymbol{\mathcal{T}}`$.
pub type CauchyTangentStiffness = TensorRank4<3, 1, 1, 1, 0>;

/// The rate tangent stiffness associated with the Cauchy stress $`\boldsymbol{\mathcal{V}}`$.
pub type CauchyRateTangentStiffness = TensorRank4<3, 1, 1, 1, 0>;

/// A list of coordinates.
pub type Coordinates<const I: usize, const W: usize> = TensorRank1List<3, I, W>;

/// A coordinate in the current configuration.
pub type CurrentCoordinate = TensorRank1<3, 1>;

/// A list of coordinates in the current configuration.
pub type CurrentCoordinates<const W: usize> = TensorRank1List<3, 1, W>;

/// A velocity in the current configuration.
pub type CurrentVelocity = TensorRank1<3, 1>;

/// The deformation gradient $`\mathbf{F}`$.
pub type DeformationGradient = TensorRank2<3, 1, 0>;

/// The deformation gradient rate $`\dot{\mathbf{F}}`$.
pub type DeformationGradientRate = TensorRank2<3, 1, 0>;

/// A list of deformation gradients.
pub type DeformationGradients<const W: usize> = TensorRank2List<3, 1, 0, W>;

/// A list of deformation gradient rates.
pub type DeformationGradientRates<const W: usize> = TensorRank2List<3, 1, 0, W>;

/// The first Piola-Kirchoff stress $`\mathbf{P}`$.
pub type FirstPiolaKirchoffStress = TensorRank2<3, 1, 0>;

/// A list of first Piola-Kirchoff stresses.
pub type FirstPiolaKirchoffStresses<const W: usize> = TensorRank2List<3, 1, 0, W>;

/// The tangent stiffness associated with the first Piola-Kirchoff stress $`\boldsymbol{\mathcal{C}}`$.
pub type FirstPiolaKirchoffTangentStiffness = TensorRank4<3, 1, 0, 1, 0>;

/// A list of first Piola-Kirchoff tangent stiffnesses.
pub type FirstPiolaKirchoffTangentStiffnesses<const W: usize> = TensorRank4List<3, 1, 0, 1, 0, W>;

/// The rate tangent stiffness associated with the first Piola-Kirchoff stress $`\boldsymbol{\mathcal{U}}`$.
pub type FirstPiolaKirchoffRateTangentStiffness = TensorRank4<3, 1, 0, 1, 0>;

/// A list of first Piola-Kirchoff rate tangent stiffnesses.
pub type FirstPiolaKirchoffRateTangentStiffnesses<const W: usize> = TensorRank4List<3, 1, 0, 1, 0, W>;

/// A force.
pub type Force = TensorRank1<3, 1>;

/// A list of forces.
pub type Forces<const W: usize> = TensorRank1List<3, 1, W>;

/// The frame spin $`\mathbf{\Omega}=\dot{\mathbf{Q}}\cdot\mathbf{Q}^T`$.
pub type FrameSpin = TensorRank2<3, 1, 1>;

/// The heat flux.
pub type HeatFlux = TensorRank1<3, 1>;

/// The left Cauchy-Green deformation $`\mathbf{B}`$.
pub type LeftCauchyGreenDeformation = TensorRank2<3, 1, 1>;

/// A coordinate in the reference configuration.
pub type ReferenceCoordinate = TensorRank1<3, 0>;

/// A list of coordinates in the reference configuration.
pub type ReferenceCoordinates<const W: usize> = TensorRank1List<3, 0, W>;

/// The right Cauchy-Green deformation $`\mathbf{C}`$.
pub type RightCauchyGreenDeformation = TensorRank2<3, 0, 0>;

/// The rotation of the current configuration $`\mathbf{Q}`$.
pub type RotationCurrentConfiguration = TensorRank2<3, 1, 1>;

/// The rate of rotation of the current configuration $`\dot{\mathbf{Q}}`$.
pub type RotationRateCurrentConfiguration = TensorRank2<3, 1, 1>;

/// The rotation of the reference configuration $`\mathbf{Q}_0`$.
pub type RotationReferenceConfiguration = TensorRank2<3, 0, 0>;

/// A scalar.
pub type Scalar = TensorRank0;

/// A list of scalars.
pub type Scalars<const W: usize> = TensorRank0List<W>;

/// The second Piola-Kirchoff stress $`\mathbf{S}`$.
pub type SecondPiolaKirchoffStress = TensorRank2<3, 0, 0>;

/// The tangent stiffness associated with the second Piola-Kirchoff stress $`\boldsymbol{\mathcal{G}}`$.
pub type SecondPiolaKirchoffTangentStiffness = TensorRank4<3, 0, 0, 1, 0>;

/// The rate tangent stiffness associated with the second Piola-Kirchoff stress $`\boldsymbol{\mathcal{W}}`$.
pub type SecondPiolaKirchoffRateTangentStiffness = TensorRank4<3, 0, 0, 1, 0>;

/// A stiffness resulting from a force.
pub type Stiffness = TensorRank2<3, 1, 1>;

/// A list of stiffnesses.
pub type Stiffnesses<const W: usize> = TensorRank2List2D<3, 1, 1, W>;

/// The temperature gradient.
pub type TemperatureGradient = TensorRank1<3, 1>;

/// A vector.
pub type Vector<const I: usize> = TensorRank1<3, I>;

/// A list of vectors.
pub type Vectors<const I: usize, const W: usize> = TensorRank1List<3, I, W>;

/// A 2D list of vectors.
pub type Vectors2D<const I: usize, const W: usize, const X: usize> = TensorRank1List2D<3, I, W, X>;
