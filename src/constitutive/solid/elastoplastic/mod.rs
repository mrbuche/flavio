//! Elastoplastic constitutive models.

#[cfg(test)]
pub mod test;

use super::
{
    *,
    elastic::Elastic,
    super::fluid::plastic::Plastic
};

/// Required methods for elastoplastic constitutive models.
pub trait Elastoplastic<'a>
where
    Self: Elastic<'a> + Plastic<'a>
{
    fn compute_elastic_deformation(&self, deformation_gradient: &DeformationGradient) -> DeformationGradientElastic
    {
        deformation_gradient * self.get_plastic_deformation_gradient().inverse()
    }
    fn get_plastic_deformation_gradient(&self) -> &DeformationGradientPlastic;
    fn set_plastic_deformation_gradient(&mut self, plastic_deformation_gradient: DeformationGradientPlastic);
}
