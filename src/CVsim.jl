module CVsim

  using Parameters
  using QuantumOptics
  using Permutations
  using LinearAlgebra

  # Parameters of the simulations
  export CVsim_Parameters
  include("parameters.jl")

  # Fast calculation of the Wigner quasiprobability
  export wignerify
  include("wignerify.jl")

  # State preparation functions
  export prepare_GKP
  export prepare_Cat
  include("prepare_GKP.jl")
  include("prepare_Cat.jl")

  # State stabilization functions
  export stabilize_GKP
  include("stabilize_GKP.jl")
  
end