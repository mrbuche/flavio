# Arruda-Boyce model

**Parameters**

- The bulk modulus $\kappa$.
- The shear modulus $\mu$.
- The number of links $N_b$.

**State variables**

- The deformation gradient $\mathbf{F}$.

**Notes**

- The nondimensional end-to-end length per link of the chains is $\gamma=\sqrt{\mathrm{tr}(\mathbf{B}^*)/3N_b}$.
- The nondimensional force is given by the inverse Langevin function as $\eta=\mathcal{L}^{-1}(\gamma)$.
- The initial values are given by $\gamma_0=\sqrt{1/3N_b}$ and $\eta_0=\mathcal{L}^{-1}(\gamma_0)$.
- The Arruda-Boyce model reduces to the [Neo-Hookean model](neo-hookean.md) when $N_b\to\infty$.

**Helmholtz free energy density**

$$
a = \frac{3\mu N_b\gamma_0}{\eta_0}\left[\gamma\eta - \gamma_0\eta_0 - \ln\left(\frac{\eta_0\sinh\eta}{\eta\sinh\eta_0}\right) \right] + \frac{\kappa}{2}\left[\frac{1}{2}\left(J^2 - 1\right) - \ln J\right]
$$

**Cauchy stress**

$$
\boldsymbol{\sigma} = \frac{\mu\gamma_0\eta}{J\gamma\eta_0}\,{\mathbf{B}^*}' + \frac{\kappa}{2}\left(J - \frac{1}{J}\right)\mathbf{1}
$$

**Tangent stiffness**

$$
\mathcal{T}_{ijkL} = \frac{\mu\gamma_0\eta}{J^{5/3}\gamma\eta_0}\left(\delta_{ik}F_{jL} + \delta_{jk}F_{iL} - \frac{2}{3}\,\delta_{ij}F_{kL}- \frac{5}{3} \, B_{ij}'F_{kL}^{-T} \right) + \frac{\kappa}{2} \left(J + \frac{1}{J}\right)\delta_{ij}F_{kL}^{-T} \\ + \frac{\mu\gamma_0\eta}{3J^{7/3}N_b\gamma^2\eta_0}\left(\frac{1}{\eta\mathcal{L}'(\eta)} - \frac{1}{\gamma}\right)B_{ij}'B_{km}'F_{mL}^{-T}
$$
