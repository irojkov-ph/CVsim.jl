"""
    stabilize_GKP( rho, P; method = "DeNeeve", kwargs...)

Stabilization of Gottesman-Kitaev-Preskill (GKP) states using various methods

Inputs:
  - rho::Operator{B,B}    Operator (density matrix) defined in the basis B 
  - P::CVsim_Parameters   Parameter structure as defined in parameters.jl
  - method::String        Method to prepare a GKP state.
                          Default: "DeNeeve"
  - kwargs                Other parameters specific to each method. 
Output:
  - A GKP state as an Operator (density matrix)

"""
function stabilize_GKP(rho::Operator, 
                       P::CVsim_Parameters;
                       method::String = "DeNeeve", kwargs...)
  if method == "DeNeeve"
    f_name = "DeNeeve_stabilization"
  elseif method == "CampagneIbarcq"
    f_name = "CampagneIbarcq_stabilization"
  elseif method == "Royer"
    f_name = "Royer_stabilization"
  elseif method == "Sivak"
    f_name = "Sivak_stabilization"
  end

  # Generate the correct f-function
  f = getfield(CVsim, Symbol(f_name))

  return f(rho, P; kwargs...)

end;

"""
    stabilize_GKP( psi, P; method = "DeNeeve", kwargs...)

Stabilization of Gottesman-Kitaev-Preskill (GKP) states using various methods
(for now only De Neeve's method is available).

Inputs:
  - psi::StateVector{B}   State vector (Ket or Bra) defined in the basis B
  - P::CVsim_Parameters   Parameter structure as defined in parameters.jl
  - method::String        Method to prepare a GKP state.
                          Default: "DeNeeve" 
  - kwargs                Other parameters specific to each method. 
Output:
  - A GKP state as an Operator (density matrix)

"""
function stabilize_GKP(psi::StateVector, 
                       P::CVsim_Parameters; 
                       method::String = "DeNeeve", kwargs...)

  return stabilize_GKP(dm(psi), P; method=method, kwargs...)
end;

"""
    DeNeeve_stabilization( rho, P; n_times = 1, ε = 2√π*0.045, α = √π, μ = 2√π*0.065, option = "", indices = 0, kraus = true)

Stabilization of Gottesman-Kitaev-Preskill (GKP) using the spin
degree of freedom of an ion as reported by De Neeve et al. in
DOI: 10.1038/s41567-021-01487-7

Inputs:
  - rho::Operator{B,B}    Operator (density matrix) defined in the 
                          position/momentum basis B or a composite basis
                          B with at least one position/momentum basis in it. 
  - P::CVsim_Parameters   Parameter structure as defined in parameters.jl
  - n_times::Int          Number of stabilization rounds
  - ε, α, μ ::Float       Three parametres as defined in the paper
  - option::String        String specifying which can contain "q" or "p" 
                          for correcting only one of the quadratures. 
                          If none, both will be stabilized. 
  - indices               Index or vector of indices of subspaces that one 
                          would like to stabilize. If index = 0, then all 
                          the relevant subspaces will be subjected to n_cycles 
                          stabilization.
  - kraus::Bool           Use Kraus operators instead of the unitaries and the
                          spin degree of freedom. Default: true.
Output:
  - A GKP state as an Operator (density matrix)

"""
function DeNeeve_stabilization(rho::Operator, P::CVsim_Parameters, index::Int;
                               n_times::Int = 1,
                               ε = 2*sqrt(π)*0.045,
                               α = sqrt(π),
                               μ = 2*sqrt(π)*0.065,
                               option = "",
                               kraus = true )

  # If "q"/"p" not specified in options add both
  option = (occursin("q",option)||occursin("p",option)) ? option : option*"qp"

  # Index of the spin subspace
  if index == 0
    spin_subsp_ind = 2
  else
    num_basis = length(basis(rho).bases)
    spin_subsp_ind = num_basis + 1

    # Place the subspace that we want to stabilize 
    # at the end of the basis, we repermute the dimension 
    # after stabilization
    perm = vcat(setdiff([1:num_basis;], [index]),[index])
    p = Permutation(perm)
    
    rho = perm != [1:num_basis;] ? permutesystems(rho, perm) : rho;
  end

  # Take the operator of only the last subspace
  op_sub(op) = index == 0 ? op : one(ptrace(basis(rho),num_basis)) ⊗ op

  # Create Kraus map
  if kraus
    cos(p::Operator) = Operator(p.basis_l,p.basis_r,
                        LinearAlgebra.isdiag(p.data) ? 
                                 LinearAlgebra.Diagonal(Base.cos.(diag(p.data))) :
                                 LinearAlgebra.cos(dense(p).data) )
    sin(p::Operator) = Operator(p.basis_l,p.basis_r,
                        LinearAlgebra.isdiag(p.data) ? 
                                 LinearAlgebra.Diagonal(Base.sin.(diag(p.data))) :
                                 LinearAlgebra.sin(dense(p).data) )
    # This is a more efficient calculation of exp(1im*A) for a matrix A
    cis(p::Operator) = Operator(p.basis_l,p.basis_r,
                        LinearAlgebra.isdiag(p.data) ?
                                LinearAlgebra.Diagonal(Base.cis.(diag(p.data))) :
                                LinearAlgebra.cis(dense(p).data) )

    K(c,q,p,pm) = op_sub( 1/sqrt(2.0) * cis(pm*c[1]*p) *
                          ( cos(c[2]*q) * cis(pm*c[3]*p) + pm * 
                            sin(c[2]*q) * cis(-pm*c[3]*p) ) )
  end;

  # De Neeve operators
  # q-quadrature stabilization cycles definition
  if occursin("q",option) && ~kraus
    epsilon_gate(u_k) = exp(dense(1im * u_k * momentum(P.b) ⊗ sigmay(P.b_spin)));
    alpha_gate(v_k)   = exp(dense(1im * v_k * position(P.b) ⊗ sigmax(P.b_spin)));
    mu_gate(w_k)      = exp(dense(1im * w_k * momentum(P.b) ⊗ sigmay(P.b_spin)));
    stabq = op_sub(mu_gate(μ)*alpha_gate(α)*epsilon_gate(ε))
  elseif occursin("q",option)
    Kq_p = sparse(K( [μ,α,ε], position(P.b), momentum(P.b), +1))
    Kq_m = sparse(K( [μ,α,ε], position(P.b), momentum(P.b), -1))
  end;

  # p-quadrature stabilization cycles definition
  if occursin("p",option) && ~kraus
    epsilon2_gate(u_k) = exp(dense(1im * u_k * position(P.b) ⊗ sigmay(P.b_spin)));
    alpha2_gate(v_k)   = exp(dense(-1im * v_k * momentum(P.b) ⊗ sigmax(P.b_spin)));
    mu2_gate(w_k)      = exp(dense(1im * w_k * position(P.b) ⊗ sigmay(P.b_spin)));
    stabp = op_sub(mu2_gate(μ)*alpha2_gate(α)*epsilon2_gate(ε));
  elseif occursin("p",option)
    Kp_p = sparse(K( [μ,α,ε], -momentum(P.b), position(P.b), +1))
    Kp_m = sparse(K( [μ,α,ε], -momentum(P.b), position(P.b), -1))
  end;

  # Define the evolution
  stabq_cycle(rho) = ~kraus ?
                     ptrace(op_sub(stabq)*(rho⊗dm(spinup(P.b_spin)))*op_sub(dagger(stabq)),spin_subsp_ind) :
                     Kq_p * rho * dagger(Kq_p) + Kq_m * rho * dagger(Kq_m);
  stabp_cycle(rho) = ~kraus ?
                     ptrace(op_sub(stabp)*(rho⊗dm(spinup(P.b_spin)))*op_sub(dagger(stabp)),spin_subsp_ind) :
                     Kp_p * rho * dagger(Kp_p) + Kp_m * rho * dagger(Kp_m);

  # Run stabilization cycles
  for i = [1:n_times;]
    if occursin("q",option)
      rho = stabq_cycle(rho)
    end;
    if occursin("p",option)
      rho = stabp_cycle(rho)
    end;
  end;

  # Repermute the dimension after stabilization
  if index != 0 && perm != [1:num_basis;]
    perm = inv(p).data    
    rho = permutesystems(rho, perm)
  end

  return rho;

end;

function DeNeeve_stabilization(rho::Operator{B,B},
                               P::CVsim_Parameters; kwargs...) where B<:PositionBasis  

  QuantumOpticsBase.check_samebases(basis(rho),P.b)
  return DeNeeve_stabilization(rho, P, 0; kwargs...);

end;

function DeNeeve_stabilization(rho::Operator{B,B}, 
                               P::CVsim_Parameters; kwargs...) where B<:MomentumBasis

  QuantumOpticsBase.check_samebases(basis(rho),P.b)
  return DeNeeve_stabilization(rho, P, 0; kwargs...);

end;

function DeNeeve_stabilization(rho::Operator{B,B}, 
                               P::CVsim_Parameters;
                               indices = 0, 
                               kwargs... ) where B<:CompositeBasis

  # If the user specifies a vector of indices, then
  # execute the function for each index sequentially
  for index in indices
    # Check if the index is well defined
    basis_length = length(basis(rho).bases)
    @assert ((0 <= index <= basis_length) && isa(index,Int)) 
            "Index must be an integer between 0 and the number of subspaces."
    
    # If the index is different from 0, check if the basis 
    # of the subspace specified by the index is the same as P.b
    if index != 0
      QuantumOpticsBase.check_samebases(basis(rho).bases[index],P.b)
      return DeNeeve_stabilization(rho, P, index; kwargs...);
    end;

    # Determine the subspaces' indices with the same basis as P.b
    where_bases = [ QuantumOpticsBase.samebases(basis(rho).bases[k],P.b) ? k : 0 for k=1:basis_length ];
    deleteat!(where_bases, where_bases .== 0);

    for i in where_bases
      rho = DeNeeve_stabilization(rho, P, i; kwargs...);
    end;
  end;

  return rho;

end;

"""
    CampagneIbarcq_stabilization( rho, P; n_times = 1, l = 2√π, Δ = 0.37, option="all", indices = 0)

Stabilization of Gottesman-Kitaev-Preskill (GKP) using an auxillary
transmon and consisting of a sharpening and trimming steps for each 
quadrature as reported by Campagne-Ibarcq et al. in
DOI: 10.1038/s41586-020-2603-3

Inputs:
  - rho::Operator{B,B}    Operator (density matrix) defined in the 
                          position/momentum basis B or a composite basis
                          B with at least one position/momentum basis in it. 
  - P::CVsim_Parameters   Parameter structure as defined in parameters.jl
  - n_times::Int          Number of stabilization rounds
  - l::Float              Lattice constant of a square lattice encoding
  - Δ::Float              Δ² - variance of Gaussian peaks in the GKP state
                          1/Δ² - variance of the Gaussian envelope
  - option::String        String specifying which steps to execute:
                              "qS"     - only sharpening of q-quadrature
                              "qST"    - sharpening-trimming of q only
                          Default: "all" which is equivallent to "qpST"
  - indices               Index or vector of indices of subspaces that one 
                          would like to stabilize. If index = 0, then all 
                          the relevant subspaces will be subjected to n_cycles 
                          stabilization.
Output:
  - A GKP state as an Operator (density matrix)

"""
function CampagneIbarcq_stabilization(rho::Operator, P::CVsim_Parameters, index::Int;
                                      n_times::Int = 1,
                                      l::Real = 2*sqrt(pi), 
                                      Δ::Real = 0.37,
                                      option::String = "all" )

  # Change the option "all" to default value "qpST" + change "st" to capital letters
  option = option=="all" ? "qpST" : option
  option = String([((c == 's') || (c == 't')) ? uppercase(c) : c for c in option])

  # Define the constants as in Royer et al. paper
  cΔ = cosh(Δ^2); sΔ = sinh(Δ^2);

  # Index of the spin subspace
  if index == 0
    spin_subsp_ind = 2
  else
    num_basis = length(basis(rho).bases)
    spin_subsp_ind = num_basis + 1

    # Place the subspace that we want to stabilize 
    # at the end of the basis, we repermute the dimension 
    # after stabilization
    perm = vcat(setdiff([1:num_basis;], [index]),[index])
    p = Permutation(perm)

    rho = perm != [1:num_basis;] ? permutesystems(rho, perm) : rho;
  end

  op_sub(op) = index == 0 ? op : one(ptrace(basis(rho),num_basis)) ⊗ op
  # Campagne-Ibarcq operators
  # q-quadrature sharpening cycles definition
  if occursin("q",option) && occursin("S",option)
    sharpenq_gate  = op_sub(exp(dense(-1im * (sΔ*l/2.0) * momentum(P.b) ⊗ sigmay(P.b_spin))));
    sharpenq_gate *= op_sub(exp(dense(-1im * (cΔ*l/2.0) * position(P.b) ⊗ sigmax(P.b_spin))));
  end;
  if occursin("q",option) && occursin("T",option)
    trimq_gate  = op_sub(exp(dense(-1im * (cΔ*l/2.0) * position(P.b) ⊗ sigmax(P.b_spin))));
    trimq_gate *= op_sub(exp(dense(-1im * (sΔ*l/2.0) * momentum(P.b) ⊗ sigmay(P.b_spin))));
  end;
  if occursin("p",option) && occursin("S",option)
    sharpenp_gate  = op_sub(exp(dense( 1im * (sΔ*l/2.0) * position(P.b) ⊗ sigmay(P.b_spin))));
    sharpenp_gate *= op_sub(exp(dense(-1im * (cΔ*l/2.0) * momentum(P.b) ⊗ sigmax(P.b_spin))));
  end;
  if occursin("p",option) && occursin("T",option)
    trimp_gate  = op_sub(exp(dense(-1im * (cΔ*l/2.0) * momentum(P.b) ⊗ sigmax(P.b_spin))));
    trimp_gate *= op_sub(exp(dense( 1im * (sΔ*l/2.0) * position(P.b) ⊗ sigmay(P.b_spin))));
  end;

  # Define the plus state of the ancillary qubit used for the stabilization + spindown(P.b_spin)
  plus_st = normalize(spinup(P.b_spin))

  # Run stabilization cycles
  for i = [1:n_times;]
    if occursin("q",option) && occursin("S",option)
      rho = sharpenq_gate * ( rho ⊗ dm(plus_st) ) * dagger(sharpenq_gate) ;
      rho = ptrace(rho, spin_subsp_ind );
    end;
    if occursin("q",option) && occursin("T",option)
      rho = trimq_gate * ( rho ⊗ dm(plus_st) ) * dagger(trimq_gate) ;
      rho = ptrace(rho, spin_subsp_ind );
    end;
    if occursin("p",option) && occursin("S",option)
      rho = sharpenp_gate * ( rho ⊗ dm(plus_st) ) * dagger(sharpenp_gate) ;
      rho = ptrace(rho, spin_subsp_ind );
    end;
    if occursin("p",option) && occursin("T",option)
      rho = trimp_gate * ( rho ⊗ dm(plus_st) ) * dagger(trimp_gate) ;
      rho = ptrace(rho, spin_subsp_ind );
    end;
  end;

  # Repermute the dimension after stabilization
  if index != 0 && perm != [1:num_basis;]
    perm = inv(p).data    
    rho = permutesystems(rho, perm)
  end

  return rho;

end;

function CampagneIbarcq_stabilization(rho::Operator{B,B},
                               P::CVsim_Parameters; kwargs...) where B<:PositionBasis  

  QuantumOpticsBase.check_samebases(basis(rho),P.b)
  return CampagneIbarcq_stabilization(rho, P, 0; kwargs...);

end;

function CampagneIbarcq_stabilization(rho::Operator{B,B}, 
                                      P::CVsim_Parameters; kwargs...) where B<:MomentumBasis

  QuantumOpticsBase.check_samebases(basis(rho),P.b)
  return CampagneIbarcq_stabilization(rho, P, 0; kwargs...);

end;

function CampagneIbarcq_stabilization(rho::Operator{B,B}, 
                                      P::CVsim_Parameters;
                                      indices = 0, 
                                      kwargs... ) where B<:CompositeBasis

  # If the user specifies a vector of indices, then
  # execute the function for each index sequentially
  for index in indices
    # Check if the index is well defined
    basis_length = length(basis(rho).bases)
    @assert ((0 <= index <= basis_length) && isa(index,Int)) 
    "Index must be an integer between 0 and the number of subspaces."

    # If the index is different from 0, check if the basis 
    # of the subspace specified by the index is the same as P.b
    if index != 0
      QuantumOpticsBase.check_samebases(basis(rho).bases[index],P.b)
      return CampagneIbarcq_stabilization(rho, P, index; kwargs...);
    end;

    # Determine the subspaces' indices with the same basis as P.b
    where_bases = [ QuantumOpticsBase.samebases(basis(rho).bases[k],P.b) ? k : 0 for k=1:basis_length ];
    deleteat!(where_bases, where_bases .== 0);

    for i in where_bases
      rho = CampagneIbarcq_stabilization(rho, P, i; kwargs...);
    end;
  end;

  return rho;

end;


"""
    Royer_stabilization( rho, P; n_times = 1, l = 2√π, Δ = 0.37, option="sBs", indices = 0)

Stabilization of Gottesman-Kitaev-Preskill (GKP) using an auxillary
qubit and either of the two methods Big-small-Big or small-Big-small
as reported by Royer et al. in Eqs. (6b) and (6c), respectively, from
DOI: 10.1038/s41586-020-2603-3

Inputs:
  - rho::Operator{B,B}    Operator (density matrix) defined in the 
                          position/momentum basis B or a composite basis
                          B with at least one position/momentum basis in it. 
  - P::CVsim_Parameters   Parameter structure as defined in parameters.jl
  - n_times::Int          Number of stabilization rounds
  - l::Float              Lattice constant of a square lattice encoding
  - Δ::Float              Δ² - variance of Gaussian peaks in the GKP state
                          1/Δ² - variance of the Gaussian envelope
  - option::String        String specifying which method to use:
                            "BsB" - Big-small-Big (work well for Δ ≤ 0.25)
                            "sBs" - small-Big-small (work better in general)
                          The option can also contain "q" or "p" for correcting
                          only one of the quadratures. If none, both will be
                          stabilized.
                          Default: "sBs"
  - indices               Index or vector of indices of subspaces that one 
                          would like to stabilize. If index = 0, then all 
                          the relevant subspaces will be subjected to n_cycles 
                          stabilization.
Output:
  - A GKP state as an Operator (density matrix)

"""
function Royer_stabilization(rho::Operator, P::CVsim_Parameters, index::Int;
                             n_times::Int = 1,
                             l::Real = 2*sqrt(pi), 
                             Δ::Real = 0.37,
                             option::String="sBs",
                             kraus::Bool = false )

  @assert (occursin("BsB",option) || occursin("sBs",option)) 
          "The option must contain at least BsB or sBs."

  # If "q"/"p" not specified in options add both
  option = (occursin("q",option)||occursin("p",option)) ? option : option*"qp"

  # Define the constants as in Royer et al. paper
  cΔ = cosh(Δ^2); sΔ = sinh(Δ^2);

  # Index of the spin subspace
  if index == 0
    spin_subsp_ind = 2
  else
    num_basis = length(basis(rho).bases)
    spin_subsp_ind = num_basis + 1

    # Place the subspace that we want to stabilize 
    # at the end of the basis, we repermute the dimension 
    # after stabilization
    perm = vcat(setdiff([1:num_basis;], [index]),[index])
    p = Permutation(perm)

    rho = perm != [1:num_basis;] ? permutesystems(rho, perm) : rho;
  end

  op_sub(op) = index == 0 ? op : one(ptrace(basis(rho),num_basis)) ⊗ op

  # Create Kraus map if the kraus option is marked true
  if kraus
    sqrt_I_sin(p::Operator,pm) = 
        Operator(p.basis_l,p.basis_r,
                 LinearAlgebra.sqrt(LinearAlgebra.I - pm*LinearAlgebra.sin(dense(p).data)))

    K(ϵ,q,p,pm) = op_sub( 0.5 * exp( 1im*l*cΔ*abs(ϵ)/8 - 1im*π/2) * sqrt_I_sin((l*cΔ*q-ϵ/2*p),pm)  +
                          0.5 * exp(-1im*l*cΔ*abs(ϵ)/8 + 1im*π/2) * sqrt_I_sin((l*cΔ*q+ϵ/2*p),pm) )
  end;

  # Royer sBs operators
  if occursin("sBs",option)
    if occursin("q",option) && ~kraus
      stabq  = exp(dense(-1im * (sΔ*l/4.0) * momentum(P.b) ⊗ sigmay(P.b_spin)));
      stabq *= exp(dense(-1im * (cΔ*l/2.0) * position(P.b) ⊗ sigmax(P.b_spin)));
      stabq *= exp(dense(-1im * (sΔ*l/4.0) * momentum(P.b) ⊗ sigmay(P.b_spin)));
    elseif occursin("q",option)
      Kq_p = sparse( K(l*sΔ, position(P.b), momentum(P.b), +1) )
      Kq_m = sparse( K(l*sΔ, position(P.b), momentum(P.b), -1) )
    end;
    if occursin("p",option) && ~kraus 
      stabp  = exp(dense(-1im * (sΔ*l/4.0) * position(P.b) ⊗ sigmay(P.b_spin)));
      stabp *= exp(dense( 1im * (cΔ*l/2.0) * momentum(P.b) ⊗ sigmax(P.b_spin)));
      stabp *= exp(dense(-1im * (sΔ*l/4.0) * position(P.b) ⊗ sigmay(P.b_spin)));
    elseif occursin("p",option)
      Kp_p = sparse( K(l*sΔ, -momentum(P.b), position(P.b), +1) )
      Kp_m = sparse( K(l*sΔ, -momentum(P.b), position(P.b), -1) )
    end;
  # Royer BsB operators
  else
    if occursin("q",option) && ~kraus
      stabq  = exp(dense(-1im * (cΔ*l/2.0) * position(P.b) ⊗ sigmax(P.b_spin)));
      stabq *= exp(dense(-1im * (sΔ*l) * momentum(P.b) ⊗ sigmay(P.b_spin)));
      stabq *= exp(dense(-1im * (cΔ*l/2.0) * position(P.b) ⊗ sigmax(P.b_spin)));
    elseif occursin("q",option)
      Kq_p = sparse( K(-4*l*sΔ, position(P.b), momentum(P.b), -1) )
      Kq_m = sparse( K(-4*l*sΔ, position(P.b), momentum(P.b), +1) )
    end;
    if occursin("p",option) && ~kraus
      stabp  = exp(dense( 1im * (cΔ*l/2.0) * momentum(P.b) ⊗ sigmax(P.b_spin)));
      stabp *= exp(dense(-1im * (sΔ*l) * position(P.b) ⊗ sigmay(P.b_spin)));
      stabp *= exp(dense( 1im * (cΔ*l/2.0) * momentum(P.b) ⊗ sigmax(P.b_spin)));
    elseif occursin("p",option)
      Kp_p = sparse( K(-4*l*sΔ, -momentum(P.b), position(P.b), -1) )
      Kp_m = sparse( K(-4*l*sΔ, -momentum(P.b), position(P.b), +1) )
    end;
  end;

  # Define the evolution
  stabq_cycle(rho) = ~kraus ?
                     ptrace(op_sub(stabq)*(rho⊗dm(spinup(P.b_spin)))*op_sub(dagger(stabq)),spin_subsp_ind) :
                     Kq_p * rho * dagger(Kq_p) + Kq_m * rho * dagger(Kq_m);
  stabp_cycle(rho) = ~kraus ?
                     ptrace(op_sub(stabp)*(rho⊗dm(spinup(P.b_spin)))*op_sub(dagger(stabp)),spin_subsp_ind) :
                     Kp_p * rho * dagger(Kp_p) + Kp_m * rho * dagger(Kp_m);

  # Run stabilization cycles
  for i = [1:n_times;]
    if occursin("q",option)
      rho = stabq_cycle(rho)
    end;
    if occursin("p",option)
      rho = stabp_cycle(rho)
    end;
  end;

  # Repermute the dimension after stabilization
  if index != 0 && perm != [1:num_basis;]
    perm = inv(p).data    
    rho = permutesystems(rho, perm)
  end

  return rho;

end;

function Royer_stabilization(rho::Operator{B,B},
                             P::CVsim_Parameters; kwargs...) where B<:PositionBasis  

  QuantumOpticsBase.check_samebases(basis(rho),P.b)
  return Royer_stabilization(rho, P, 0; kwargs...);

end;

function Royer_stabilization(rho::Operator{B,B}, 
                             P::CVsim_Parameters; kwargs...) where B<:MomentumBasis

  QuantumOpticsBase.check_samebases(basis(rho),P.b)
  return Royer_stabilization(rho, P, 0; kwargs...);

end;

function Royer_stabilization(rho::Operator{B,B}, 
                             P::CVsim_Parameters;
                             indices = 0, 
                             kwargs... ) where B<:CompositeBasis

  # If the user specifies a vector of indices, then
  # execute the function for each index sequentially
  for index in indices
    # Check if the index is well defined
    basis_length = length(basis(rho).bases)
    @assert ((0 <= index <= basis_length) && isa(index,Int)) 
    "Index must be an integer between 0 and the number of subspaces."

    # If the index is different from 0, check if the basis 
    # of the subspace specified by the index is the same as P.b
    if index != 0
      QuantumOpticsBase.check_samebases(basis(rho).bases[index],P.b)
      return Royer_stabilization(rho, P, index; kwargs...);
    end;

    # Determine the subspaces' indices with the same basis as P.b
    where_bases = [ QuantumOpticsBase.samebases(basis(rho).bases[k],P.b) ? k : 0 for k=1:basis_length ];
    deleteat!(where_bases, where_bases .== 0);

    for i in where_bases
      rho = Royer_stabilization(rho, P, i; kwargs...);
    end;
  end;

  return rho;

end;




"""
    Sivak_stabilization( rho, P; n_times = 1, l = 2√π, Δ = 0.37, option="sBs", indices = 0)

Stabilization of Gottesman-Kitaev-Preskill (GKP) code using the same 
methods as in Royer et al. (2020) and De Neeve et al. (2022) but to 
the exception that here the stabilizer measurement outcomes are not 
discarded but recorded. The simulation is thus a state vector simulation
which outputs the bosonic state as a state vector and the measurement outcome.
DOI: 10.1038/s41586-023-05782-6

Inputs:
  - rho::Ket{B}           State vector defined in the position/momentum 
                          basis B or a composite basis B with at least 
                          one position/momentum basis in it. 
  - P::CVsim_Parameters   Parameter structure as defined in parameters.jl
  - n_times::Int          Number of stabilization rounds
  - ε, α, μ ::Float       Three parametres as defined DeNeeve_stabilization
  - option::String        String specifying which can contain "q" or "p" 
                          for correcting only one of the quadratures. 
                          If none, both will be stabilized. 
  - indices               Index or vector of indices of subspaces that one 
                          would like to stabilize. If index = 0, then all 
                          the relevant subspaces will be subjected to n_cycles 
                          stabilization.
  - kraus::Bool           Use Kraus operators instead of the unitaries and the
                          spin degree of freedom. Default: true.
Output:
  - A GKP state as a state vector (Ket)
  - Outcomes of the stabilizer measurement (Bool)

"""
function Sivak_stabilization(psi::Ket, 
                             P::CVsim.CVsim_Parameters, 
                             index::Int;
                             n_times::Int = 1,
                             ε = 2*sqrt(π)*0.045,
                             α = sqrt(π),
                             μ = 2*sqrt(π)*0.065,
                             option = "",
                             kraus = true
                            )

  # If "q"/"p" not specified in options add both
  option = (occursin("q",option)||occursin("p",option)) ? option : option*"qp"

  # Defining the syndrome vector (initially empty)
  syndromes = zeros(0)  

  # Index of the spin subspace
  if index == 0
    spin_subsp_ind = 2
  else
    num_basis = length(basis(psi).bases)
    spin_subsp_ind = num_basis + 1

    # Place the subspace that we want to stabilize 
    # at the end of the basis, we repermute the 
    # dimension after stabilization
    perm = vcat(setdiff([1:num_basis;], [index]),[index])
    p = Permutation(perm)
    
    psi = perm != [1:num_basis;] ? permutesystems(psi, perm) : psi;
  end

  # Take the operator of only the last subspace
  op_sub(op) = index == 0 ? op : one(ptrace(basis(psi),num_basis)) ⊗ op

  # Create Kraus map
  if kraus
    cos(p::Operator) = Operator(p.basis_l,p.basis_r,
                        LinearAlgebra.isdiag(p.data) ? 
                                 LinearAlgebra.Diagonal(Base.cos.(diag(p.data))) :
                                 LinearAlgebra.cos(dense(p).data) )
    sin(p::Operator) = Operator(p.basis_l,p.basis_r,
                        LinearAlgebra.isdiag(p.data) ? 
                                 LinearAlgebra.Diagonal(Base.sin.(diag(p.data))) :
                                 LinearAlgebra.sin(dense(p).data) )
    # This is a more efficient calculation of exp(1im*A) for a matrix A
    cis(p::Operator) = Operator(p.basis_l,p.basis_r,
                        LinearAlgebra.isdiag(p.data) ?
                                LinearAlgebra.Diagonal(Base.cis.(diag(p.data))) :
                                LinearAlgebra.cis(dense(p).data) )

    Kpm(c,q,p,pm) = 1/sqrt(2.0) * cis(pm*c[1]*p) *
                   ( cos(c[2]*q) * cis(pm*c[3]*p) + pm * 
                     sin(c[2]*q) * cis(-pm*c[3]*p) )

    K(c,q,p,pm) = op_sub( 1/sqrt(2.0) * (Kpm(c,q,p,+1) + pm*Kpm(c,q,p,-1)) )

  end;

  # De Neeve operators
  # q-quadrature stabilization cycles definition
  if occursin("q",option) && ~kraus
    epsilon_gate(u_k) = exp(dense(1im * u_k * momentum(P.b) ⊗ sigmay(P.b_spin)));
    alpha_gate(v_k)   = exp(dense(1im * v_k * position(P.b) ⊗ sigmax(P.b_spin)));
    mu_gate(w_k)      = exp(dense(1im * w_k * momentum(P.b) ⊗ sigmay(P.b_spin)));
    stabq = op_sub(mu_gate(μ)*alpha_gate(α)*epsilon_gate(ε))
  elseif occursin("q",option)
    Kq_p = sparse(K( [μ,α,ε], position(P.b), momentum(P.b), +1))
    Kq_m = sparse(K( [μ,α,ε], position(P.b), momentum(P.b), -1))
  end;

  # p-quadrature stabilization cycles definition
  if occursin("p",option) && ~kraus
    epsilon2_gate(u_k) = exp(dense(1im * u_k * position(P.b) ⊗ sigmay(P.b_spin)));
    alpha2_gate(v_k)   = exp(dense(-1im * v_k * momentum(P.b) ⊗ sigmax(P.b_spin)));
    mu2_gate(w_k)      = exp(dense(1im * w_k * position(P.b) ⊗ sigmay(P.b_spin)));
    stabp = op_sub(mu2_gate(μ)*alpha2_gate(α)*epsilon2_gate(ε));
  elseif occursin("p",option)
    Kp_p = sparse(K( [μ,α,ε], -momentum(P.b), position(P.b), +1))
    Kp_m = sparse(K( [μ,α,ε], -momentum(P.b), position(P.b), -1))
  end;

  # Define the evolution
  stabq_cycle(psi) = 
      ~kraus ?
              measurement_process(op_sub(stabq)*(psi⊗spinup(P.b_spin)),
                                  [ op_sub(identityoperator(P.b)⊗dm(spinup(P.b_spin))),
                                    op_sub(identityoperator(P.b)⊗dm(spindown(P.b_spin))) ]) :
              measurement_process(psi,[Kq_p,Kq_m]);
  stabp_cycle(psi) = 
      ~kraus ?
              measurement_process(op_sub(stabp)*(psi⊗spinup(P.b_spin)),
                                  [ op_sub(identityoperator(P.b)⊗dm(spinup(P.b_spin))),
                                    op_sub(identityoperator(P.b)⊗dm(spindown(P.b_spin))) ]) :
              measurement_process(psi,[Kp_p,Kp_m]);

  # Run stabilization cycles
  for i = [1:n_times;]
    if occursin("q",option)
      psi, s = stabq_cycle(psi)
      append!(syndromes,s)
    end;
    if occursin("p",option)
      psi, s = stabp_cycle(psi)
      append!(syndromes,s)
    end;
  end;

  # Repermute the dimension after stabilization
  if index != 0 && perm != [1:num_basis;]
    perm = inv(p).data    
    psi = permutesystems(psi, perm)
  end

  return (psi,syndromes);

end;

function Sivak_stabilization(psi::Ket{B}, 
                             P::CVsim.CVsim_Parameters;
                             indices = 0, 
                             kwargs... ) where B<:CompositeBasis

  # Defining the syndrome vector (initially empty)
  syndromes = zeros(0)

  # If the user specifies a vector of indices, then
  # execute the function for each index sequentially
  for index in indices
    # Check if the index is well defined
    basis_length = length(basis(psi).bases)
    @assert ((0 <= index <= basis_length) && isa(index,Int)) 
    "Index must be an integer between 0 and the number of subspaces."

    # If the index is different from 0, check if the basis 
    # of the subspace specified by the index is the same as P.b
    if index != 0
      QuantumOpticsBase.check_samebases(basis(psi).bases[index],P.b)
      return DeNeeve_stabilization(psi, P, index; kwargs...);
    end;

    # Determine the subspaces' indices with the same basis as P.b
    where_bases = [ QuantumOpticsBase.samebases(basis(psi).bases[k],P.b) ? k : 0 for k=1:basis_length ];
    deleteat!(where_bases, where_bases .== 0);

    for i in where_bases
      psi, s = DeNeeve_stabilization(psi, P, i; kwargs...);
      append!(syndromes,s)
    end;

    syndromes = reshape(syndromes,div(length(syndromes),length(where_bases)),length(where_bases))

  end;

  return (psi,syndromes);

end;

function Sivak_stabilization(psi::Ket, 
                             P::CVsim.CVsim_Parameters;
                             indices = 0,
                             kwargs... )

  return DeNeeve_stabilization(psi, P, 0; kwargs...);
end;

# Helper function used in Sivak_stabilization
function measurement_process(psi::Ket, operators)

  # calculate the probabilities
  proba = ones(length(operators))
  for (i,op) in enumerate(operators)
    tmp = op*psi
    proba[i] = abs(dagger(tmp)*(tmp))
  end;

  # check the normalization
  proba_new = sum(proba) == 1.0 ? proba : proba ./ sum(proba)

  # sample a single operator from the list of 
  # operators according to the probabilities
  Op = sample(operators,Weights(proba_new))
  syndrome = findall(x -> x == Op, operators)
  
  # apply the operator
  psi = normalize(Op*psi)

  return (psi,syndrome)
end;
