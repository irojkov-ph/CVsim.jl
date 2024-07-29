# [Wigner quasiprobability of a single CV system state](@id wigner_single)

The [Wigner quasiprobability](https://en.wikipedia.org/wiki/Wigner_quasiprobability_distribution) representation of a single continuous-variable (CV) system's state ``\rho`` is defined as

```math
  W(x,p)=\frac{1}{\pi\hbar}\int_{-\infty}^{+\infty} \,
  \exp\left(\frac{2i}{\hbar}\,p\,y\right) \,
  \langle x-y \vert \rho \vert x+y\rangle \,
  \mathrm{d}y \,.
```

It is evident that this representation is defined by the Fourier transform of the antidiagonals of the density matrix ``\rho`` expressed in the position basis. Since `CVsim.jl` relies on position or momentum representations of CV system states, we can obtain their Wigner quasiprobability through a fast Fourier transform of the columns of their "skewed density matrix". The latter corresponds to the original position density matrix but skewed along one direction,

```math
  \begin{pmatrix}
  \cdot & \cdot & \cdot & \cdot \\
  \cdot & \cdot & \cdot & \cdot \\
  \cdot & \cdot & \cdot & \cdot \\
  \cdot & \cdot & \cdot & \cdot 
  \end{pmatrix} 
  \mapsto
  \begin{pmatrix}
  \cdot & \cdot & \cdot & \cdot & 0 & 0 & 0 \\
  0 & \cdot & \cdot & \cdot & \cdot & 0 & 0\\
  0 & 0 & \cdot & \cdot & \cdot & \cdot & 0\\
  0 & 0 & 0 & \cdot & \cdot & \cdot & \cdot 
\end{pmatrix} 
```

Note that the first, third, fifth, etc. columns will contain the desired antidiagonal elements. The even columns will contain in-between antidiagonals and hold mostly redundant information. They could be included but that requires extra computation due to how they straddle across the diagonal. Selecting the odd columns and performing their FFT gives us the desired Wigner quasi-probability as a square matrix of the same dimension as ``\rho``.

This calculation has been first implemented in Python by the Wikipedia user [Nanite](https://en.wikipedia.org/wiki/User:Nanite) in their ["Long-time evolution of a mixed state ρ in an anharmonic  potential well."](https://en.wikipedia.org/wiki/File:Hamiltonian_flow_quantum.webm) simulation. We adapted it in Julia here.

```@docs
wignerify(::Matrix{ComplexF64})
```

In `CVsim.jl`, we overload this function to enable seamless integration with `QuantumOptics.jl`, specifically `Operator` objects (representing density matrices) and `StateVector` objects (encompassing `Ket` and `Bra` states), all defined in a position basis.

```@docs
wignerify(rho::Operator{B,B}) where B<:PositionBasis
wignerify(psi::StateVector{B}) where B<:PositionBasis
```

