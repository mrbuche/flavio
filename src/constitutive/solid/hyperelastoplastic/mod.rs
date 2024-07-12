//! Hyperelastoplastic constitutive models.
//!
//! ---
//!
//! Hyperelastoplastic constitutive models are defined by a Helmholtz free energy density function and a plastic work density function.
//!
//! ```math
//! \mathbf{P}:\dot{\mathbf{F}} - \dot{a}(\mathbf{F}^\mathrm{e}) - \dot{w}(\mathbf{F}^\mathrm{p}) \geq 0
//! ```
//! The two functions are separate functions of the elastic and plastic deformation gradients, respectively, which are related to the applied deformation gradient by the Bilby-Kr&ouml;ner-Lee decomposition. The plastic work density function is rarely specified directly.
//!
//! ```math
//! \mathbf{F} = \mathbf{F}^\mathrm{e}\cdot\mathbf{F}^\mathrm{p}
//! ```
//! The plastic deformation gradient is selected as the internal state variable, and plastic deformation is assumed to be incompressible. Extremizing the sum of the Helmholtz free energy density and plastic work density functions with respect to the plastic deformation gradient then shows that the plastic Cauchy stress is equal to the elastic Mandel stress.
//!
//! ```math
//! \boldsymbol{\sigma}^\mathrm{p} = \mathbf{M}^\mathrm{e} = J\,{\mathbf{F}^\mathrm{e}}^T\mathbf{F}^T\cdot\boldsymbol{\sigma}^\mathrm{e}\cdot{\mathbf{F}^\mathrm{e}}^{-T}
//! ```
//! Incompressible plasticity implies a deviatoric plastic velocity gradient and therefore deviatoric plastic stretching and spin tensors. Accordingly, the plastic stress is assumed to be deviatoric. The plastic flow is assumed irrotational, so the plastic spin tensor is zero. Since the elastic and plastic deformations can evolve independently, the second law of thermodynamics can be satisfied separately. The plastic stretching rate becomes proportional to the deviatoric Mandel stress, often written in terms of a scalar plastic flow rate.
//!
//! ```math
//! \mathbf{D}^\mathrm{p} = \dot{\gamma}^\mathrm{p}\,\frac{{\mathbf{M}^\mathrm{e}}'}{|{\mathbf{M}^\mathrm{e}}'|}
//! ```
//! The stress becomes hyperelastic, where the Cauchy stress can be shown to equal the elastic Cauchy stress.
//!
//!```math
//! \mathbf{P} = \frac{\partial a}{\partial\mathbf{F}} = \mathbf{P}^\mathrm{e}\cdot{\mathbf{F}^\mathrm{p}}^{-T}
//! ```
//! Consequently, the tangent stiffness associated with the first Piola-Kirchoff stress is symmetric for hyperelastoplastic models.
//!
//! ```math
//! \mathcal{C}_{iJkL} = \mathcal{C}_{kLiJ}
//! ```

#[cfg(test)]
pub mod test;

mod saint_venant_kirchoff;

pub use saint_venant_kirchoff::SaintVenantKirchoff;

use super::
{
    *,
    elastic::Elastic,
    elastoplastic::Elastoplastic,
    hyperelastic::Hyperelastic,
    super::fluid::plastic::Plastic,
};

/// Required methods for hyperelastoplastic constitutive models.
pub trait Hyperelastoplastic<'a>
where
    Self: Elastoplastic<'a> + Hyperelastic<'a>
{}
