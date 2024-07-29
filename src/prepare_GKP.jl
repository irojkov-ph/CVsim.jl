using Theta
using Distributions

"""
    prepare_GKP( P; method = "default", kwargs...)

Preparation of Gottesman-Kitaev-Preskill (GKP) states using various methods.

Inputs:
  - `P::CVsim_Parameters`   : Parameter structure as defined in parameters.jl
  - `method::String`        : Method to prepare a GKP state. Currently only
                              "default", "Hastrup", "Albert", "Mensen" methods  
                              are implemented. See their definitions bellow.
  - `kwargs`                : Other parameters specific to each method. 
Output:
  - A GKP state as a `Ket` state or as an `Operator` depending on the method (density matrix)

"""
function prepare_GKP(P::CVsim_Parameters;
                     method::String = "default", kwargs...)
  if method == "default"
      f_name = "default_GKP_preparation"
  elseif method == "Hastrup"
    f_name = "hastrup_preparation"
  elseif method == "Albert"
    f_name = "albert_preparation"
  elseif method == "Mensen"
    f_name = "mensen_preparation"
  else
    throw(ArgumentError("Method not implemented."))
  end

  # Generate the correct f-function
  f = getfield(CVsim, Symbol(f_name))

  return f(P; kwargs...)

end;

"""
    hastrup_preparation( P; n_cycles = 2, κ_squeez = .37, σ_ground = 1.0)

Preparation of the approximate GKP state using an additional two-level
quantum system and measurement-free preparation protocol as proposed by
Hastrup et al. in DOI: 10.1038/s41534-020-00353-3

Warning: The returned state has a spin component.
Inputs:
  - `P::CVsim_Parameters`  : Parameter structure as defined in parameters.jl
  - `n_cycles::Int`        : Number of cycles of preparation as defined in 
                             the paper. Suported numbers are 2, 3, 4.
  - `κ_squeez::Float64`    : Squeezing of the initial squeezed vacuum state
  - `σ_ground::Float64`    : Standard deviation of the vacuum state (gaussian state)
Output:
  - |0⟩L GKP Ket ⊗ spin Ket

"""
function hastrup_preparation(P::CVsim_Parameters;
                             n_cycles::Int = 2,
                             κ_squeez::Float64 = .37,
                             σ_ground::Float64 = 1.0)

  @assert 2<=n_cycles<=4 "In Hastrup preparation, n_cycles should be 2, 3 or 4."

  # Squeezing operator
  r = -log(κ_squeez);
  Sq = exp(dense( 1im * r/2 * (position(P.b)*momentum(P.b) + momentum(P.b)*position(P.b)) ));
  
  # Create the initial squeezed vacuum state needed for the protocol
  x0 = 0.0; p0 = 0.0
  ψ_init = Sq*gaussianstate(P.b, x0, p0, σ_ground) ⊗ spinup(P.b_spin);
  
  # Hastrup operators
  preparation_gate(u_k)   = exp(dense(1im * u_k * position(P.b) ⊗ sigmay(P.b_spin) ));
  displacement_gate(v_k)  = exp(dense(1im * v_k * momentum(P.b) ⊗ sigmax(P.b_spin) ));
  disentangling_gate(w_k) = exp(dense(1im * w_k * position(P.b) ⊗ sigmay(P.b_spin) ));

  prep_cycle(u_k,v_k,w_k) = disentangling_gate(w_k)*displacement_gate(v_k)*preparation_gate(u_k);

  # Constants for each prepatation cycle (u_k taken from paper) 
  u = [ [0.0, 0.045,   0.0, 0.0],
        [0.0, 0.053, 0.033, 0.0],
        [0.0, 0.038, 0.027, 0.015] ]

  function v(k,N)
      if k==1
          return -sqrt(pi) * 2.0^(N-1)
      else
          return sqrt(pi) * 2.0^(N-k)
      end
  end;

  function w(k,N)
      if k<N
          return -sqrt(pi)/4 * 2.0^(-(N-k))
      else
          return sqrt(pi)/4
      end
  end;

  # Three prepatation cycles
  GKP_0 = ψ_init;

  for k=1:n_cycles
      u_k = u[n_cycles-1][k+1];
      v_k = v(k,n_cycles);
      w_k = w(k,n_cycles);
      GKP_0 = prep_cycle(u_k,v_k,w_k)*GKP_0
  end

  return GKP_0

end;

"""
    albert_preparation( P; Δ = 0.37, n_peaks = 51, GKP_1 = false)  

Preparation of the approximate GKP state by summing a certain number of
displaced squeezed states as reported by Albert et al. in Eq.(7.7a) from
DOI: 10.1103/PhysRevA.97.032346 

Inputs:
  - `P::CVsim_Parameters`   : Parameter structure as defined in parameters.jl
  - `Δ::Float64`            : Δ² - variance of Gaussian peaks in the GKP state,
                              and 1/Δ² - variance of the Gaussian envelope
  - `n_peaks::Int`          : Number of peaks to initialize
  - `GKP_1::Bool`           : if set to true |1⟩L will be initialized otherwise |0⟩L
Output:
  - |GKP⟩ Ket vector

"""
function albert_preparation(P::CVsim_Parameters;
                            Δ::Float64 = 0.37,
                            n_peaks::Int = 51,
                            GKP_1::Bool = false)

    displace(α) = exp(-1im*imag(α)*real(α)) * 
                  LazyProduct( exp(1im*√2.0 * dense( imag(α)*position(P.b))), P.Tpm,
                               exp(1im*√2.0 * dense(-real(α)*momentum(P.b_momentum))), P.Tmp )

    squeeze(Δ)  = exp(1im * Δ/2 * dense( 1im * identityoperator(P.b) + 2 * momentum(P.b)*position(P.b) ))

    ground = gaussianstate(P.b, 0.0, 0.0, 1)

    GKP_0 = 0*gaussianstate(P.b, 0.0, 0.0, 1)

    if mod(n_peaks,2)==0
        n_peaks += 1
    end

    n = - Int((n_peaks-1)/2)
    for i in [0:n_peaks;]
        GKP_0 = GKP_0 + exp(-(π/2) * Δ^2 * (2*n+GKP_1)^2 ) * displace(√(π/2)*(2*n+GKP_1)) * squeeze(-log(Δ)) * ground
        n=n+1
    end

    return normalize(GKP_0)

end;

"""
    mensen_preparation( P; μ = 0.0, Δ = 0.37, κ = 0.37, period = 2√π, logic = [[0],[0]]) 

Preparation of the approximate GKP state using the Jacobi θ function as
as reported by Mensen et al. in Eq.(3.34) from
DOI: 10.1103/PhysRevA.104.022408

To evaluate the Jacobi θ function we use the package 
[Theta.jl](https://github.com/chualynn/Theta.jl) from
DOI: 10.2140/jsag.2021.11.41

Inputs:
  - `P::CVsim_Parameters`        : Parameter structure as defined in parameters.jl
  - `μ::Float64`                 : Center of the big gaussian envelope
  - `Δ::Float64`                 : Δ² - variance of Gaussian peaks in the GKP state
  - `κ::Float64`                 : 1/κ² - variance of the Gaussian envelope
  - `period::Float64`            : Period of the θ function, 
    * set 2√π for |0⟩L and |1⟩L states
    * set √π for |+⟩L and |-⟩L states
  - `logic::Vector{Vector{Int}}` : Vector of two booleans which sets the order of the peaks
    * set [[0],[0]] for |0⟩L and |+⟩L
    * set [[0],[1]] for |1⟩L
    * set [[1],[0]] for |-⟩L
Output:
  - |GKP⟩ Ket vector

"""
function mensen_preparation(P::CVsim_Parameters;
                            μ::Float64 = 0.0, 
                            Δ::Float64 = 0.37,
                            κ::Float64 = 0.37,
                            period::Float64 = 2*sqrt(pi),
                            logic::Vector{Vector{Int}}=[[0],[0]]) 
  
  normalization = sqrt(4π*Δ/κ)
  envelope = pdf.(Normal(μ,κ^(-1)), P.qvec)
  R = RiemannMatrix([2π*1im*Δ^2.0 /(period^2)])
  f(x)= 1/(sqrt(period)) * theta([x / period], R, char=logic)
  wavefunction = normalization * envelope .* f.(P.qvec)

  GKP = Ket(P.b,wavefunction);

  return GKP

end;


"""
    default_GKP_preparation( P; μ = 0.0, Δ = 0.37, κ = 0.37, period = 2√π, parity = 0) 

Preparation of the approximate GKP state using a sum of gaussians with an outer envelope

Inputs:
  - `P::CVsim_Parameters` : Parameter structure as defined in parameters.jl
  - `μ::Float64`          : Center of the big gaussian envelope
  - `Δ::Float64`          : Δ² - variance of Gaussian peaks in the GKP state
  - `κ::Float64`          : 1/κ² - variance of the Gaussian envelope
  - `period::Float64`     : Period of the θ function, 
    * set 2√π for |0⟩L and |1⟩L states
    * set √π for |+⟩L and |-⟩L states
  - `parity::Bool`        : Boolean that sets the parity of the state
    * set 0 for |0⟩L and |+⟩L
    * set 1 for |1⟩L and |-⟩L
  - `n_peaks::Int`        : Half of the total number of Gaussian peaks that
                            that compose the state
Output:
  - |GKP⟩ Ket vector

"""
function default_GKP_preparation(P::CVsim_Parameters;
                                 μ::Float64 = 0.0, 
                                 Δ::Float64 = 0.37,
                                 κ::Float64 = 0.37,
                                 period::Float64 = 2*sqrt(pi),
                                 parity::Bool = false,
                                 n_peaks::Int = 20
                                ) 
  
  # Determine if one wants to create a minus state
  minus_state = ( parity && (period == sqrt(pi)) )
  parity_term(ind) = (-1)^(ind*minus_state)
  
  one_state = ( parity && (period == 2*sqrt(pi)) )
  envelope(ind) = exp(- κ^2/2 * ( period*ind - μ + period/2*Int(one_state))^2 )
  peaks(x,ind) = exp(- (2*Δ^2)^(-1) * (x - period*ind - period/2*Int(one_state))^2)
  
  wavefunction = sum([ parity_term(ind).*envelope(ind).*peaks.(P.qvec,ind) 
                       for ind=[-n_peaks+round(μ):1:n_peaks+round(μ);] ] )

  GKP = Ket(P.b,convert(Vector{ComplexF64},wavefunction));

  return normalize(GKP)

end;