# [General structure of state stabilization scripts](@id stabilize_general)

`CVsim.jl` was designed for the simulation of continuous variable systems in quadrature space with a particular emphasize on the simulation of bosonic encodings. The goal of the state stabilization scripts is to implement various stabilization protocols for bosonic codes. Some of the stabilization schemes have been experimentally implemented while others are only theoretical. We focussed so far only on measurement-free [stabilization protocols for GKP states](@ref stabilize_GKP). Moreover, the implemented schemes are not continuous in time and correspond in practice to a sequence of pulses/gates. It would be interesting to implement additional dissipative discrete/continuous stabilization schemes as well as measurement-based error-correction protocols for GKP codes and other rotation/translation-symmetric bosonic encodings such as Cat and binomial codes. Multi-mode codes are also in principle implementable in `CVsim.jl`. For a detailed overview of various existing bosonic codes we recommend to have a look at the [Error Correction Zoo](https://errorcorrectionzoo.org/).

The general structure of the state preparation scripts is the following:
The filename of the script must be `stabilize_NAME.jl` with `NAME` denoting the states that are meant to be stabilized with this script. The two main method definitions must then be

```julia
  stabilize_NAME(rho::Operator, P::CVsim_Parameters; method::String = "default", kwargs...)
  stabilize_NAME(psi::StateVector, P::CVsim_Parameters; method::String = "default", kwargs...)
```

Unlike the [state stabilization scripts](@ref prepare_general), here the main method must take exactly two arguments: First, the input state of type `StateVector` (i.e. `Ket` or `Bra`) or `Operator` (i.e. a density matrix); second, `P` corresponding to the parameters of a given continuous variable simulation (see [CVsim Parameters](@ref parameters)). The optional variable `method` then allows the user to select the stabilization scheme/protocol that they desire and `kwargs` correspond to some additional parameters that this scheme requires.

!!! info
    Feel free to contribute in adding novel stabilization protocols following this general structure. If you think there could be a better way to implement this in Julia feel free to propose it (we are not experts in Julia). In both cases, please [open a pull request](https://github.com/irojkov-ph/CVsim.jl/pulls) with your proposed additions/modifications.

