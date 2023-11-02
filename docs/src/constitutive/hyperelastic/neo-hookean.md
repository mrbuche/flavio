# Neo-Hookean model

**Parameters**

- The bulk modulus $\kappa$.
- The shear modulus $\mu$.

**State variables**

- The deformation gradient $\mathbf{F}$.

**Helmholtz free energy density**

$$
a = \frac{\mu}{2}\left[\mathrm{tr}(\mathbf{B}^*) - 3\right] + \frac{\kappa}{2}\left[\frac{1}{2}\left(J^2 - 1\right) - \ln J\right]
$$

**Cauchy stress**

$$
\boldsymbol{\sigma} = \frac{\mu}{J}\,{\mathbf{B}^*}' + \frac{\kappa}{2}\left(J - \frac{1}{J}\right)\mathbf{1}
$$

**Tangent stiffness**

$$
\mathcal{T}_{ijkL} = \frac{\mu}{J^{5/3}}\left(\delta_{ik}F_{jL} + \delta_{jk}F_{iL} - \frac{2}{3}\,\delta_{ij}F_{kL}- \frac{5}{3} \, B_{ij}'F_{kL}^{-T} \right) + \frac{\kappa}{2} \left(J + \frac{1}{J}\right)\delta_{ij}F_{kL}^{-T}
$$
