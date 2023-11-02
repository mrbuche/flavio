# Mooney-Rivlin model

**Parameters**

- The bulk modulus $\kappa$.
- The shear modulus $\mu$.
- The extra modulus $\mu_m$.

**State variables**

- The deformation gradient $\mathbf{F}$.

**Notes**

- The Mooney-Rivlin model reduces to the [Neo-Hookean model](neo-hookean.md) when $\mu_m\to 0$.

**Helmholtz free energy density**

$$
a = \frac{\mu - \mu_m}{2}\left[\mathrm{tr}(\mathbf{B}^* ) - 3\right] + \frac{\mu_m}{2}\left[I_2(\mathbf{B}^*) - 3\right] + \frac{\kappa}{2}\left[\frac{1}{2}\left(J^2 - 1\right) - \ln J\right]
$$

**Cauchy stress**

$$
\boldsymbol{\sigma} = \frac{\mu - \mu_m}{J}\, {\mathbf{B}^* }' - \frac{\mu_m}{J}\left(\mathbf{B}^{* -1}\right)' + \frac{\kappa}{2}\left(J - \frac{1}{J}\right)\mathbf{1}
$$

**Tangent stiffness**

$$
\mathcal{T}_{ijkL} = \frac{\mu-\mu_m}{J^{5/3}}\left(\delta_{ik}F_{jL} + \delta_{jk}F_{iL} - \frac{2}{3}\,\delta_{ij}F_{kL}- \frac{5}{3} \, B_{ij}'F_{kL}^{-T} \right) + \frac{\kappa}{2} \left(J + \frac{1}{J}\right)\delta_{ij}F_{kL}^{-T} \\ - \frac{\mu_m}{J}\left[ \frac{2}{3}\,B_{ij}^{* -1}F_{kL}^{-T} - B_{ik}^{* -1}F_{jL}^{-T} - B_{ik}^{* -1}F_{iL}^{-T} + \frac{2}{3}\,\delta_{ij}\left(B_{km}^{* -1}\right)'F_{mL}^{-T} - \left(B_{ij}^{* -1}\right)'F_{kL}^{-T} \right]
$$
