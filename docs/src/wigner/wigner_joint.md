# [Wigner quasiprobability of joint CV systems state](@id wigner_joint)

The Wigner quasiprobability representation of a state ``\rho`` of a ``N``-mode system (i.e. many continuous-variable systems) is defined as[^1][^2]

```math
  W(\mathbf{x})=\frac{1}{(2\pi)^{2N}}\int_{\mathbb{R}^{2N}} \,
  \exp\!\left(-i\,\mathbf{x^\mathrm{T}\Omega\xi}\right) \,
  \mathrm{Tr}\!\left[ \rho\,\exp\!\left(i\mathbf{\hat{x}^\mathrm{T}\Omega\xi}\right) \right] \,
  \mathrm{d}^{2N}\mathbf{\xi} \,.
```

where ``\mathbf{x}`` and ``\mathbf{\xi}`` are vectors comprising real values of the ``2N`` quadratures, while ``\mathbf{\hat{x}}`` is a vector encompassing the quadrature operators themselves. Here, ``\mathbf{\Omega}`` corresponds to the symplectic form. The exponential term inside the trace can be rewritten in terms of displacement operators on the ``2N`` modes. Thus, similar method as for the single mode presented in [Single system state](@ref wigner_single) but has not been implemented in `CVsim.jl`, yet. Such quasiprobability distribution could be helpful for the diagnostic of entanglement between multiple mode as it has been demonstrated in the context of logical two-qubit gate between bosonic codewords[^3].


[^1]: [K. E. Cahill and R. J. Glauber, Phys. Rev. **177**, 1882 (1969)](https://journals.aps.org/pr/abstract/10.1103/PhysRev.177.1882)
[^2]: [C. Weedbrook _et al._ Rev. Mod. Phys. **84**, 621 (2012)](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.84.621)
[^3]: [C. Wang _et al._ Science **352**, 6289, 1087-1091 (2016)](https://doi.org/10.1126/science.aaf2941)
