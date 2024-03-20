#[cfg(test)]
mod test;

use super::*;

const G: usize = 3;
const M: usize = 2;
const N: usize = 6;
const O: usize = 6;
const P: usize = 4;
const Q: usize = 3;

pub struct Triangle<C>
{
    constitutive_models: [C; G],
    projected_gradient_vectors: ProjectedGradientVectors<G, N>,
    reference_normals: ReferenceNormals<P>,
    scaled_composite_jacobians: Scalars<G>
}

// need to calculate and store reference normals for each subtriangle during initialization