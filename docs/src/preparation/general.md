# [General structure of state preparation scripts](@id prepare_general)

`CVsim.jl` was designed for the simulation of continuous variable systems in quadrature space with a particular emphasize on the simulation of bosonic encodings. The two primary classes of bosonic codes are translation-symmetric codes and rotation-symmetric codes[^1], exemplified by Gottesman-Kitaev-Preskill (GKP)[^2] and Schrödinger's Cat[^3][^4] codes, respectively.

The goal of the state preparation scripts is to implement various state preparation protocols for these codes. Some of the state preparation schemes have been experimentally implemented while others are only theoretical (meaning that they are in general not accessible in experiments). So far, mostly [preparation scripts for GKP states](@ref prepare_GKP) have been written. [Cat state preparation](@ref prepare_Cat) includes only the default preparation. It would be interesting to implement additional schemes for those codes as well as other fancy continuous-variable states such as binomial codewords, cubic phase state and other non-Gaussian states.

The general structure of the state preparation scripts is the following:
The filename of the script must be `prepare_NAME.jl` with `NAME` denoting the states that are meant to be prepared with this script. The main method definition must then be

```julia
  prepare_NAME(P::CVsim_Parameters; method::String = "default", kwargs...)
```

That means the main method must take exactly one parameter `P` corresponding to the parameters of a given continuous variable simulation (see [CVsim Parameters](@ref parameters)). The optional variable `method` then allows the user to select the state preparation scheme/protocol that they desire and `kwargs` correspond to some additional parameters that this scheme requires.

!!! info
    Feel free to contribute in adding novel schemes following this general structure. If you think there could be a better way to implement this in Julia feel free to propose it (we are not experts in Julia). In both cases, please [open a pull request](https://github.com/irojkov-ph/CVsim.jl/pulls) with your proposed additions/modifications.

[^1]: [A. L. Grimsmo, J. Combes, and B. Q. Baragiola, Phys. Rev. X **10**, 011058 (2020)](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.011058)
[^2]: [D. Gottesman, A. Kitaev, and J. Preskill, Phys. Rev. A **64**, 012310 (2001)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.64.012310)
[^3]: [P. T. Cochrane, G. J. Milburn, and W. J. Munro, Phys. Rev. A **59**, 2631 (1999)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.59.2631)
[^4]: [M. Mirrahimi _et al._, New J. Phys. **16**, 045014 (2014)](https://iopscience.iop.org/article/10.1088/1367-2630/16/4/045014)