# Hyperelastic models

- [Arruda-Boyce model](hyperelastic/arruda-boyce.md)
- [Gent model](hyperelastic/gent.md)
- [Mooney-Rivlin model](hyperelastic/mooney-rivlin.md)
- [Neo-Hookean model](hyperelastic/neo-hookean.md)
- [Yeoh model](hyperelastic/yeoh.md)

## Helmholtz free energy density

A hyperelastic constitutive model can be defined by specifying the Helmholtz free energy density as some function of the deformation gradient

$$
a = a(\mathbf{F})
$$

## Cauchy stress

The Cauchy stress is calculated from the Helmholtz free energy density as

$$
\boldsymbol{\sigma} = \frac{1}{J}\frac{\partial a}{\partial\mathbf{F}}\cdot\mathbf{F}^T
$$

## Tangent stiffness

The tangent stiffness is calculated from the Cauchy stress as

$$
\boldsymbol{\mathcal{T}} = \frac{\partial\boldsymbol{\sigma}}{\partial\mathbf{F}}
$$
