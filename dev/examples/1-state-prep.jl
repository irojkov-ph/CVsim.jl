# # State preparation and visualization 
#
# This example shows the basic functionality of the `CVsim.jl` package.
# 
# We start by including the necessary libraries. `CVsim.jl` can be seen 
# as an extensions of the [`QuantumOptics.jl`](https://qojulia.org/) package
# which provides the basic functionality for quantum optics simulations in
# Julia. `CVsim.jl` extends this functionality to continuous-variable (CV)
# system simulations in the position/momentum basis. We then include both 
# packages.
# 

using QuantumOpticsBase
using CVsim

# !!! note 
#     Here, we include `QuantumOpticsBase.jl` as we only require the bare
#     functionalities of the `QuantumOptics.jl` package. If one needs more
#     advanced features, one can include the full package `using QuantumOptics`.
#
# The first step on a CV simulation is to define its parameters. We will use 
# 100 points in the interval ``[-6\sqrt{\pi},6\sqrt{\pi}]`` for the
# quadrature variables. The accuracy of the simulation can be adjusted by
# changing the number of points and the size of the interval.

## CVsim parameters
P = CVsim_Parameters(n_points = 100, qvec_min = -6*√π, qvec_max= 6*√π);

# The main quadrature basis accessed using `P.b` is the position one, but
# it can be adjusted using the parameter `b_prefered` (see
# [Parameters](@ref parameters)). Using this preferred quadrature basis
# we can define the various operators such as the quadrature operators 
# `q` and `p`.

q = position(P.b);

# Using a quadrature basis as opposed to the Fock basis allows us to define
# some useful operators such as the displacement operator `D(α)` in a more 
# convenient way as some of them are diagonal in the position quadrature, e.g.

using LinearAlgebra
displace_q(α::Real) = Operator(P.b, Diagonal(exp.(1im*α*diag(q.data))));

# The momentum operator in the position basis is defined using
p = momentum(P.b);

# However, ``p`` is related to ``q`` by a Fourier transform, thus we can 
# switch between the the momentum and position bases using the transformation
# operators `P.Tmp` and `P.Tpm`. In the momentum basis, the momentum 
# operator is diagonal, thus allowing the user to improve the efficiency
# of the simulation when working with momentum operators.
#
# We can now initialize some CV states such as Cat and GKP states. For this, 
# we only require the functions `prepare_Cat` and `prepare_GKP`.

GKP = prepare_GKP(P; κ = 0.35, Δ = 0.35);
Cat = prepare_Cat(P; α = 2.5,  θ = π/4);

# We would now like to vizualize the position distribution of those states.
# For this, we include the plotting package [`Plots.jl`](https://docs.juliaplots.org/stable/)
# as well as the [`LaTeXStrings.jl`](https://github.com/JuliaStrings/LaTeXStrings.jl) 
# package for LaTeX rendering of the labels.

using Plots
using LaTeXStrings

# We can easily access the position distribution of the states using the
# squared absolute value of the diagonal of the density matrix (or simply
# the squared absolute value of the state vector) when defined in the 
# position basis.

p1 = plot(P.qvec,real.(abs2.(Cat.data)),title="Cat",linewidth=2);
p2 = plot(P.qvec,real.(abs2.(GKP.data)),title="GKP",linewidth=2);
plot(p1,p2,layout=(1,2),size=(800,300),thickness_scaling=1.5,margin=2Plots.mm,
    xlabel=L"q",ylabel=L"\mathbb{P}(q)",legend=false,framestyle=:box,xlims=(-8,8))

# The momentum quadrature distribution can be obtained similarly by 
# transforming the state to the momentum basis.

p1 = plot(P.pvec,real.(abs2.((P.Tmp*Cat).data)),title="Cat",linewidth=2);
p2 = plot(P.pvec,real.(abs2.((P.Tmp*GKP).data)),title="GKP",linewidth=2);
plot(p1,p2,layout=(1,2),size=(800,300),thickness_scaling=1.5,margin=2Plots.mm,
    xlabel=L"p",ylabel=L"\mathbb{P}(p)",legend=false,framestyle=:box,xlims=(-8,8))

# Finally, we can also visualize the Wigner quasiprobability of the 
# states by first calculating it using the `wignerify` function (see
# [Wigner quasiprobability](@ref wigner_single)) and then plotting it
# using the `heatmap` function from `Plots.jl`.

Wig = wignerify(Cat);  
w1  = heatmap(P.qvec,P.pvec,Wig,legend=false,c=cgrad([:red, :white, :blue]),
              clim=(-maximum(Wig),maximum(Wig)),levels=200,framestyle=:box,
              aspect_ratio=1,title="Cat");
yaxis!(L"p", (-4,4)); xaxis!(L"q", (-4,4));

ticks  = [ k*√π for k = -4:2:4 ];
labels = [ string(round(k))*L"\sqrt{\pi}" for k = -4:2:4 ];
Wig = wignerify(GKP);  
w2  = heatmap(P.qvec,P.pvec,Wig,legend=false,c=cgrad([:red, :white, :blue]),
              clim=(-maximum(Wig),maximum(Wig)),levels=200,framestyle=:box,
              aspect_ratio=1,title="GKP");
yaxis!(L"p", (ticks[1],ticks[end])); yticks!(ticks,labels);
xaxis!(L"q", (ticks[1],ticks[end])); xticks!(ticks,labels);

plot(w1,w2,layout=(1,2),size=(800,400),thickness_scaling=1.5,margin=2Plots.mm)

# The resolution of the Wigner function can be adjusted by changing the
# number of points and the size of the interval in the `CVsim_Parameters`. 