"""
    CVsim_Parameters( ; qvec_min = -8*√π, qvec_max = 8*√π, n_points = 300, 
                        b_prefered = "position", spin_number = 1//2, qudit_levels = 3 )

Structure storing essential parameters that define a given simulation

Inputs (optional):
  - `qvec_min::Float64`          : Minimum position in phase space.
  - `qvec_max::Float64`          : Maximum position in phase space.
  - `n_points::Int`              : Number of points equidistantly distributed between `qvec_min` and `qvec_max`.
  - `b_prefered::String`         :  Preferred basis for state representation ("position" or "momentum").
  - `spin_number::Rational{Int}` :  Spin number for a spin auxillary system.
  - `qudit_levels::Int`          : Number for qudit levels for a qudit auxillary system.
Output:
  - A `CVsim_Parameters` structure containing these parameters as well as 
    the essential bases (position, momentum, spin, qudit).
"""
@with_kw mutable struct CVsim_Parameters
  
  #######################################################
  #--- Parameters related to the continuous variable ---#
  #######################################################

  # Number atoms, levels (per atom) and Fock states (common motional mode) to consider
  qvec_min::Float64 = -8*√π  
  qvec_max::Float64 =  8*√π
  n_points::Int     = 300 

  # Position and momentum basis
  b_position = PositionBasis(qvec_min, qvec_max, n_points)
  b_momentum = MomentumBasis(b_position)
  b_prefered::String = "position"
  b = b_prefered=="position" ? b_position : b_momentum

  # Sample points (usefull for the Wigner plots later)
  qvec = samplepoints(b_position)
  pvec = samplepoints(b_momentum) # p == √(2π) * q (because conjugate variables)

  # Transformation position ⇐⇒ momentum basis
  Tmp = transform(b_momentum, b_position)
  Tpm = dagger(Tmp);

  ##########################################################
  #--- Parameters related to the spin degree of freedom ---#
  ##########################################################

  # Spin number
  spin_number::Rational{Int} = 1//2

  # Spin basis
  b_spin = SpinBasis(spin_number)

  ###########################################################
  #--- Parameters related to the qudit degree of freedom ---#
  ###########################################################

  # Number of qudit levels
  qudit_levels::Int = 3

  # Qudit basis
  b_qudit = NLevelBasis(qudit_levels)

end