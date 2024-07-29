"""
    prepare_Cat( P; method = "default", kwargs...)

Preparation of Schrödinger's Cat states using various methods.

Inputs:
  - `P::CVsim_Parameters`   : Parameter structure as defined in parameters.jl
  - `method::String`        : Method to prepare a Cat state. Currently only "default" 
                              method is implemented. See its definitions bellow.
  - `kwargs`                : Other parameters specific to each method. 
Output:
  - A Cat state as a `Ket` state or as an `Operator` depending on the method (density matrix)

"""
function prepare_Cat(P::CVsim_Parameters;
                     method::String = "default", kwargs...)
  if method == "default"
    f_name = "default_2Cat_preparation"
  else
    throw(ArgumentError("Method not implemented."))
  end

  # Generate the correct f-function
  f = getfield(CVsim, Symbol(f_name))

  return f(P; kwargs...)

end;


"""
    default_2Cat_preparation( P; α = 1.0+0.0im, logical_basis = false, θ = π/2, ϕ = 0.0) 

Preparation of two-component Cat states as a superposition of two gaussians 
equally distant from the origin of phase space. The state is defined as
```
    |Cat⟩ ∝ cos(θ)|0⟩L + exp(iϕ)sin(θ)|1⟩L
```
Inputs:
  - `P::CVsim_Parameters` : Parameter structure as defined in parameters.jl
  - `α::Number`           : Amplitude of the coherent states (can be `Complex`)
  - `logical_basis::Bool` : Logical basis in which the state is prepared 
    * set false for defining |0⟩L (|1⟩L) as right (left) coherent states
    * set true  for defining |0⟩L (|1⟩L) as even (odd) superposition of coherent states
  - `θ::Float64`          : sin(θ) defines the proportion of |0⟩L and |1⟩L states
  - `ϕ::Float64`          : Relative phase between |0⟩L and |1⟩L
Output:
  - |Cat⟩ Ket vector (normalized)

"""
function default_2Cat_preparation(P::CVsim_Parameters;
                                  α::Number=1.0+0.0im, 
                                  logical_basis::Bool=false,
                                  θ::Float64 = π/2,
                                  ϕ::Float64 = 0.0
                                ) 

  zeroL = gaussianstate(P.b,real(α),imag(α),1.0)
  oneL = gaussianstate(P.b,-real(α),-imag(α),1.0)
  
  if logical_basis
    zeroL, oneL = normalize(zeroL + oneL), normalize(zeroL - oneL)
  end 
  
  return normalize(normalize(cos(θ)*zeroL + exp(1im*ϕ)*sin(θ)*oneL))
end;