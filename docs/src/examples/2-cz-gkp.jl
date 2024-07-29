# # CZ gate between two GKP states 
#
# This is a more advanced example compared to the previous one as 
# it involves now two bosonic modes which encode two logical qubit.
# As previously we start by including the necessary libraries.

using QuantumOpticsBase
using CVsim
using LinearAlgebra
using Plots
using LaTeXStrings

# We then define the parameters of our CV simulation. We will use 
# 300 points in the interval ``[-12\sqrt{\pi},12\sqrt{\pi}]`` for the
# quadrature variables. We also define the quadrature operators for 
# the two modes as well as some placeholder functions to prepare the
# states where the last two prepare the maximally entangled 
# (a.k.a. Bell) states.

## CVsim parameters
P = CVsim_Parameters(n_points = 300, qvec_min = -12*√π, qvec_max= 12*√π)
q = position(P.b)

## State preparation functions
GKP0(κ,Δ) = prepare_GKP(P; κ = κ, Δ = Δ)  
GKP1(κ,Δ) = prepare_GKP(P; κ = κ, Δ = Δ, parity = true)
GKPp(κ,Δ) = prepare_GKP(P; κ = κ, Δ = Δ, period = √π)
GKPm(κ,Δ) = prepare_GKP(P; κ = κ, Δ = Δ, period = √π, parity = true)

Φ(κ,Δ,parity) = 1/√2 * ( GKP0(κ,Δ) ⊗ GKP0(κ,Δ) + (-1)^parity * GKP1(κ,Δ) ⊗ GKP1(κ,Δ) )
Ψ(κ,Δ,parity) = 1/√2 * ( GKP0(κ,Δ) ⊗ GKP1(κ,Δ) + (-1)^parity * GKP1(κ,Δ) ⊗ GKP0(κ,Δ) );

# As in the previous example we also define a placeholder function
# for ploting the Wigner function of the states using `Plots.jl`.

## Function which automates Wigner plots
function plot_wigner(state)
    ticks  = [ k*√π for k = -6:2:6 ];
    labels = [ string(round(k))*L"\sqrt{\pi}" for k = -6:2:6 ];
    Wig = wignerify(state);
  
    heatmap(P.qvec, P.pvec, Wig, legend=false,
            c=cgrad([:red, :white, :blue]),
            clim=(-maximum(Wig),maximum(Wig)),
            levels=200, framestyle = :box, aspect_ratio=1)
    yaxis!("p", (ticks[1],ticks[end]))
    yticks!(ticks,labels)
    xaxis!("q", (ticks[1],ticks[end]))
    xticks!(ticks,labels)
  end;

# We would like to implement a logical ``CZ`` gate which is a two-qubit 
# gate that flips the sign of the second qubit if the first qubit is in 
# the state ``\vert1\rangle``. For GKP codes it is know that two-qubit
# gates correspond to specific quadrature-quadrature operators. 
# In particular, controlled phase gate is defined as
# ```math
#       CZ = \exp\left(-i q_1 q_2 \right)  
# ```
# where ``q_1`` and ``q_2`` are the quadrature operators of
# the two modes.
#
# Since `CVsim.jl` represents quantum states in position/momentum space
# we can directly implement the above operator since it is diagonal in
# the position basis

CZ(θ::Number) = SparseOperator(P.b ⊗ P.b, Diagonal(exp.(-1im*θ*diag((q ⊗ q).data))));

# We allow the operator to be parametrized by an angle ``\theta`` which
# would correspond to the application time of the quadrature-quadrature
# coupling. The actual CZ gate is implemented when ``\theta=1``.
#
# ## Truth table of the CZ gate
#
# We now would like to construct the truth table of the logical CZ gate
# for the GKP states. For this we evaluate the overlap between an input
# state ``\vert\Psi_\mathrm{in}\rangle`` and an output state 
# ``\vert\Psi_\mathrm{out}\rangle`` after applying the CZ gate, i.e.
# ```math
#     O = \langle\Psi_\mathrm{in}\vert CZ \vert\Psi_\mathrm{out}\rangle\,.
# ```
# Note that here ``O`` does not exactly correspond to the fidelity
# which would be the squared norm of this overlap. We do this calculation 
# for all the logical computational basis states as input and output
# states, that we then arrange in a truth table. 

comput_basis(κ,Δ) = [GKP0(κ,Δ)⊗GKP0(κ,Δ), GKP0(κ,Δ)⊗GKP1(κ,Δ), GKP1(κ,Δ)⊗GKP0(κ,Δ), GKP1(κ,Δ)⊗GKP1(κ,Δ)];

# We construct such truth tables for GKP states with parameters
# ``\kappa_1=\Delta_1=0.35`` and ``\kappa_2=\Delta_2=0.25``.

κ1=Δ1=0.35
κ2=Δ2=0.25

basis1 = comput_basis(κ1,Δ1)
basis2 = comput_basis(κ2,Δ2)

TT1 = basis1'.*([CZ(1.0)].*basis1)
TT2 = basis2'.*([CZ(1.0)].*basis2);

# Plotting the truth tables as heatmap with color bar limited to [-1,1],
# we get the following results

gradient=cgrad([:black,:blue,:white,:red,:black], [0.15, 0.5, 0.85]);

h1=heatmap(real(TT1),c=gradient,clim=(-1,1),framestyle=:box,aspect_ratio=1,title=L"\kappa_1=\Delta_1=0.35");
annotate!([(j, i, text(round(real(TT1)[i,j],digits=3), 8,"Computer Modern",:gray)) for i in 1:4 for j in 1:4]);

h2=heatmap(real(TT2),c=gradient,clim=(-1,1),framestyle=:box,aspect_ratio=1,title=L"\kappa_2=\Delta_2=0.25");
annotate!([(j, i, text(round(real(TT2)[i,j],digits=3), 8,"Computer Modern",:gray)) for i in 1:4 for j in 1:4]);

plot(h1,h2,layout=(1,2),size=(900,400),thickness_scaling=1.5, margin=2Plots.mm,
     xlims=(0.5,4.5),xticks = (1:4,["|00⟩","|01⟩","|10⟩","|11⟩"]),
     ylims=(0.5,4.5),yticks = (1:4,["|00⟩","|01⟩","|10⟩","|11⟩"]))

# The truth tables show that the CZ gate works as expected for the GKP
# states. The gate flips the sign of the second qubit if the first qubit
# is in the state ``\vert1\rangle``. The gate is also symmetric with
# respect to the input and output states.
# 
#
# ## Wigner quasiprobability before and after
#
# However, the CZ gate is not perfect as the overlap is does not exactly
# equal to 1. This is due to the fact that the CZ gate is not a perfect
# two-qubit gate for finite-energy GKP states. To show why we can plot
# the single mode Wigner quasiprobability of the input and output states 
# of the CZ gate.

Ψ_in  = GKP0(0.3,0.3)⊗GKP0(0.3,0.3)
Ψ_out = CZ(1.0)*Ψ_in

w1=plot_wigner(ptrace(Ψ_in,2),);
w2=plot_wigner(ptrace(Ψ_out,2));
plot(w1,w2,layout=(1,2),size=(900,400),
     thickness_scaling=1.5,margin=2Plots.mm,
     title=["Input state" "Output state"])

# We observe that the Wigner quasiprobability of the output state is slightly
# distorted compared to the input state one. This finite-energy effect then 
# affects the overlap between the input and output states of the CZ gate.
# 
# ## Physical fidelity of the CZ gate
#
# Moreover, we see from the truth tables that this overlap differs for different
# energy envelopes of the GKP states, thus the physical fidelity of the gate.
# We study this effect in more detail in the next example where we evaluate 
# the fidelity of the CZ gate for different energy envelopes of the GKP states.

Δvec = [0.15:0.01:0.4;]
states = [GKP0(Δ,Δ)⊗GKP0(Δ,Δ) for Δ in Δvec]
fidelity = abs2.(diag(states'.*([CZ(1.0)].*states)))

scatter(Δvec,fidelity,size=(700,300),legend=false,framestyle=:box,
        thickness_scaling=1.5,margin=2Plots.mm,yaxis=((0.6,0.7),0.6:0.02:0.7),
        xlabel=L"\Delta",ylabel="Fidelity")

# The plot shows that the fidelity of the CZ gate converges to the finite value
# of 0.64 for the GKP states with ``\Delta\ll1``. The reasons for these effects
# and the solutions to counteract them are discussed in more detail in Ref.[^1].
#
# [^1]: [I. Rojkov, P. M. Röggla, M. Wagener, M. Fontboté-Schmidt, S. Welte, J. Home, F. Reiter, arXiv:quantum-ph/2305.05262 (2023)](https://doi.org/10.48550/arXiv.2305.05262).
