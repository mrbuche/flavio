//! Hyperelastoplastic constitutive models.
//!
//! ---
//!
//! Hyperelastoplastic constitutive models are defined by a Helmholtz free energy density function and a plastic work density function.
//!
//! ```math
//! \mathbf{P}:\dot{\mathbf{F}} - \dot{a}(\mathbf{F}_\mathrm{e}) - \dot{w}(\mathbf{F}_\mathrm{p}) \geq 0
//! ```
//! The two functions are separate functions of the elastic and plastic deformation gradients, respectively, which are related to the applied deformation gradient by the Bilby-Kr&ouml;ner-Lee decomposition. The plastic work density function is rarely specified directly.
//!
//! ```math
//! \mathbf{F} = \mathbf{F}_\mathrm{e}\cdot\mathbf{F}_\mathrm{p}
//! ```
//! The plastic deformation gradient is selected as the internal state variable, and plastic deformation is assumed to be incompressible. Extremizing the sum of the Helmholtz free energy density and plastic work density functions with respect to the plastic deformation gradient then shows that the plastic Cauchy stress is equal to the elastic Mandel stress.
//!
//! ```math
//! \boldsymbol{\sigma}_\mathrm{p} = \mathbf{M}_\mathrm{e} = J\,\mathbf{F}_\mathrm{e}^T\cdot\boldsymbol{\sigma}_\mathrm{e}\cdot\mathbf{F}_\mathrm{e}^{-T}
//! ```
//! Incompressible plasticity implies a deviatoric plastic velocity gradient and therefore deviatoric plastic stretching and spin tensors. Accordingly, the plastic stress is assumed to be deviatoric. The plastic flow is assumed irrotational, so the plastic spin tensor is zero.
//!
//! ```math
//! ???
//! ```
//! parts can evolve separately? satisfy plastic part for rate proportional to Mandel deviator, often written "in direction" with scalar plastic rate.
//!
//! ```math
//! \mathbf{D}_\mathrm{p} = \dot{\gamma}_\mathrm{p}\,\frac{\mathbf{M}_\mathrm{e}'}{|\mathbf{M}_\mathrm{e}'|}
//! ```
//! and then the elastic part by taking hyperelasticity (and Cauchy stresses equal)
//!
//!```math
//! \mathbf{P} = \frac{\partial a}{\partial\mathbf{F}} = \mathbf{P}_\mathrm{e}\cdot\mathbf{F}_\mathrm{p}^{-T}
//! ```
//! Consequently, the tangent stiffness associated with the first Piola-Kirchoff stress is symmetric for hyperelastoplastic models.
//!
//! ```math
//! \mathcal{C}_{iJkL} = \mathcal{C}_{kLiJ}
//! ```

#[cfg(test)]
pub mod test;

use super::
{
    *,
    elastoplastic::Elastoplastic,
};

/// Required methods for hyperelastoplastic constitutive models.
pub trait Hyperelastoplastic<'a>
where
    Self: Elastoplastic<'a>
{
    /// Calculates and returns the Helmholtz free energy density.
    ///
    /// ```math
    /// a = a(\mathbf{F})
    /// ```
    fn calculate_helmholtz_free_energy_density(&self, deformation_gradient: &DeformationGradient) -> Scalar;
}
