# Gent model

**Parameters**

- The bulk modulus $\kappa$.
- The shear modulus $\mu$.
- The extensibility $J_m$.

**State variables**

- The deformation gradient $\mathbf{F}$.

**Notes**

- The Gent model reduces to the [Neo-Hookean model](neo-hookean.md) when $J_m\to\infty$.

**Helmholtz free energy density**

$$
a = -\frac{\mu J_m}{2}\,\ln\left[1 - \frac{\mathrm{tr}(\mathbf{B}^* ) - 3}{J_m}\right] + \frac{\kappa}{2}\left[\frac{1}{2}\left(J^2 - 1\right) - \ln J\right]
$$

**Cauchy stress**

$$
\boldsymbol{\sigma} = \frac{J^{-1}\mu J_m {\mathbf{B}^* }'}{J_m - \mathrm{tr}(\mathbf{B}^* ) + 3} + \frac{\kappa}{2}\left(J - \frac{1}{J}\right)\mathbf{1}
$$

**Tangent stiffness**

$$
\mathcal{T}_{ijkL} = \frac{J^{-5/3}\mu J_m}{J_m - \mathrm{tr}(\mathbf{B}^* ) + 3}\Bigg[ \delta_{ik}F_{jL} + \delta_{jk}F_{iL} - \frac{2}{3}\,\delta_{ij}F_{kL} + \frac{2{B_{ij}^* }' F_{kL}}{J_m - \mathrm{tr}(\mathbf{B}^* ) + 3} \\ - \left(\frac{5}{3} + \frac{2}{3}\frac{\mathrm{tr}(\mathbf{B}^* )}{J_m - \mathrm{tr}(\mathbf{B}^* ) + 3}\right) J^{2/3} {B_{ij}^* }' F_{kL}^{-T} \Bigg] + \frac{\kappa}{2} \left(J + \frac{1}{J}\right)\delta_{ij}F_{kL}^{-T}
$$
