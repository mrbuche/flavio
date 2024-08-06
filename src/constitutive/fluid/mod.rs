//! Fluid constitutive models.

pub mod viscous;

use crate::mechanics::Scalar;
use super::Constitutive;

/// Required methods for fluid constitutive models.
pub trait Fluid<'a>
where
    Self: Constitutive<'a>
{}
