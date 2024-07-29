using FastTransforms

"""
    wignerify(mat)

Compute the Wigner transform of position-basis density matrix which is 
a Fourier transform of the anti-diagonals of the matrix (cf. Wikipedia).

Pad the matrix along one dimension, then skew the matrix. This converts
the antidiagonals to straight columns.

```     
       |‾‾‾‾|          |‾‾‾‾000|
       |    |    -->   |0    00|
       |    |          |00    0|
       |____|          |000____|
```

Note that the first, third, fifth, etc. columns will contain the desired
antidiagonal elements. The even columns will contain in-between antidiagonals.
Here we are interested only in the former ones, i.e. the skew matrix is
a square matrix of the same dimensions as the original density matrix in the 
position basis.

This calculation has been first implemented in Python by the Wikipedia user 
"Nanite" in their "Long-time evolution of a mixed state ρ in an anharmonic 
potential well." simulation and adapted in Julia here.
src: https://en.wikipedia.org/wiki/File:Hamiltonian_flow_quantum.webm

Inputs:
  - `mat::Matrix{ComplexF64}` : Matrix representing an operator (density 
                                matrix) in a position basis.
Output:
  - `Matrix{Float64}`         : Wigner transform of the matrix 
                                (same shape as `mat`).

"""
function wignerify(mat::Matrix{ComplexF64})

  # Some constants
  wig_Np = length(mat[:,1])
  wig_pshift = wig_Np ÷ 2

  # Extract anti-diagonals from mat and organizing them in columns
  mat_skew = zeros(ComplexF64,(wig_Np, wig_Np))
  for i=1:wig_Np
    mat_skew[i, (i÷2+1):(i÷2)+(wig_pshift)] .= mat[i,mod(i+1,2)+1:2:end]
  end

  # Prefactor exp(2π𝒾k) where k is the position index (i.e., x_k=x[k]) 
  shift_amt  = ([1:wig_Np;]) .* ones(wig_Np)'
  mat_skew .*= exp.((2im*pi) .* shift_amt)
  
  # Preparation of the nonuniform fast Fourier transform
  # (see FastTransforms.jl for more details)
  scale_fft = 2.0*[0:wig_Np-1;]
  fft_op = plan_nufft1(scale_fft,eps())
  
  # Fourier transform along the columns
  wigner = mapslices(x->fft_op*x, mat_skew; dims=1)

  # Execute a linear phase shift which is equivalent to a pixel shift on
  # mat_skew (but the even antidiagonals get a half-odd pixel shift)
  shift_amt = ([0:wig_Np-1;] .- (wig_pshift)) .* [0:2:2*wig_Np-2;]'
  wigner .*= exp.((2im*pi / wig_Np) .* shift_amt)
  
  # the residual imaginary part is only numerical error
  wigner = real.(wigner)
  
  # Return normalized.
  return (1/pi) .* wigner
end;

"""
    wignerify(rho)
    wignerify(rho,indices)

Compute the Wigner quasiprobability distribution for the given operator rho.

If its basis `B` is a single position basis it returns a single distribution.

If instead `B` is a composite basis then it will calculate the Wigner
quasiprobability distribution for all the subsystems which have a position basis.

Inputs:
  - `rho::Operator{B,B}`        : Operator (density matrix) defined in 
                                  the basis `B`
  - `indices::Vector{Int}`      : Indices of the subspaces for which you
                                  would like a Wigner quasiprobability
Output:
  * if `B == PositionBasis`:
    - `Matrix{Float64}`         : Wigner quasiprobability of dimension `n = length(B)`
  * if `B == CompositeBasis` or if indices are specified:
    - `Vector{Matrix{Float64}}` : Vector with Wigner quasiprobability of
                                  dimension `n_k = length(B_k)` where `B_k`
                                  is the position basis of a subspace of `B`
  
"""
wignerify(rho::Operator{B,B}) where B<:PositionBasis = wignerify(rho.data)

wignerify(rho::Operator{B,B}, index::Int) where B<:CompositeBasis = wignerify(rho,[index])

function wignerify(rho::Operator{B,B}, indices::Vector{Int}) where B<:CompositeBasis

  # Number of basis in the composite basis
  nb_bases = length(basis(rho).bases)
  @assert all(1 .<= indices .<= nb_bases) "Indices must be between 1 and the number of subspaces."

  # Create indices to trace out
  traceout_idx = [1:nb_bases;]
  deleteat!(traceout_idx,indices)

  # Trace out the useless subsystems
  subrho = ptrace(rho,traceout_idx)

  # Calculating the Wigner distributions for the remaining subsystems
  wignerify(subrho)
end;

function wignerify(rho::Operator{B,B}) where B<:CompositeBasis

  # Number of basis in the composite basis
  nb_bases = length(basis(rho).bases)
  
  # Find subspaces that are have position bases
  pos_bases = [ any(fieldnames(typeof(basis(rho).bases[i])) .== :xmin)
                for i = 1:nb_bases ]
  pos_bases_idx = [1:nb_bases;]
  pos_bases_idx = pos_bases_idx[pos_bases]

  @assert length(pos_bases_idx)>0 "No subsystem with a PositionBasis."

  # Create the vector of Wigner distributions
  wig_vec = Vector{Matrix{Float64}}(undef,nb_bases)

  # Fill in the vector at the right indices
  f(i) = wignerify(rho,[i])
  wig_vec[pos_bases_idx] = f.(pos_bases_idx)

  return wig_vec
end;

"""
    wignerify(psi)
    wignerify(psi,indices)

Compute the Wigner quasiprobability distribution for the given state vector 
(i.e. Ket or Bra). The function just transform the state vector to a density
matrix and calculate the distribution.

Inputs:
- `psi::StateVector{B}`       : State vector (Ket or Bra) defined in 
                                the basis `B`
- `indices::Vector{Int}`      : Indices of the subspaces for which you
                                would like a Wigner quasiprobability
Output:
  * if `B == PositionBasis`:
    - `Matrix{Float64}`         : Wigner quasiprobability of dimension `n = length(B)`
  * if `B == CompositeBasis` or if indices are specified:
    - `Vector{Matrix{Float64}}` : Vector with k Wigner quasiprobability of
                                  dimension `n_k = length(B_k)` where `B_k`
                                  is the position basis of a subspace of `B`
  
"""
wignerify(psi::StateVector{B}) where B<:PositionBasis  = wignerify(dm(psi).data)
wignerify(psi::StateVector{B}, index::Int) where B<:CompositeBasis = wignerify(dm(psi),[index])
wignerify(psi::StateVector{B}, indices::Vector{Int}) where B<:CompositeBasis = wignerify(dm(psi),indices)
wignerify(psi::StateVector{B}) where B<:CompositeBasis = wignerify(dm(psi))
