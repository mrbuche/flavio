# Yeoh model

**Parameters**

- The bulk modulus $\kappa$.
- The shear modulus $\mu=\mu_1$.
- The extra moduli $\mu_n$ for $n=2\ldots N$.

**State variables**

- The deformation gradient $\mathbf{F}$.

**Notes**

- The Yeoh model reduces to the [Neo-Hookean model](neo-hookean.md) when $\mu_n\to 0$ for $n=2\ldots N$.

**Helmholtz free energy density**

$$
a = \sum_{n=1}^N \frac{\mu_n}{2}\left[\mathrm{tr}(\mathbf{B}^* ) - 3\right]^n + \frac{\kappa}{2}\left[\frac{1}{2}\left(J^2 - 1\right) - \ln J\right]
$$

**Cauchy stress**

$$
\boldsymbol{\sigma} = \sum_{n=1}^N \frac{n\mu_n}{J}\left[\mathrm{tr}(\mathbf{B}^* ) - 3\right]^{n-1}\,{\mathbf{B}^*}' + \frac{\kappa}{2}\left(J - \frac{1}{J}\right)\mathbf{1}
$$

**Tangent stiffness**

$$
\mathcal{T}_{ijkL} = \sum_{n=1}^N \frac{n\mu_n}{J^{5/3}}\left[\mathrm{tr}(\mathbf{B}^* ) - 3\right]^{n-1}\left(\delta_{ik}F_{jL} + \delta_{jk}F_{iL} - \frac{2}{3}\,\delta_{ij}F_{kL}- \frac{5}{3} \, B_{ij}'F_{kL}^{-T} \right) \\ + \sum_{n=2}^N \frac{2n(n-1)\mu_n}{J^{7/3}}\left[\mathrm{tr}(\mathbf{B}^* ) - 3\right]^{n-2}B_{ij}'B_{km}'F_{mL}^{-T} + \frac{\kappa}{2} \left(J + \frac{1}{J}\right)\delta_{ij}F_{kL}^{-T}
$$
