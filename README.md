# CVsim.jl: Continuous variable quantum system simulation package 

| **Documentation**         | **Build Status**                      |
|:------------------------- |:------------------------------------- |
| [![][docs-img]][docs-url] | [![][gh-actions-img]][gh-actions-url] |

[docs-img]: https://img.shields.io/badge/docs-latest%20release-blue.svg
[docs-url]: https://irojkov-ph.github.io/CVsim.jl/

[gh-actions-img]: https://github.com/irojkov-ph/CVsim.jl/actions/workflows/CI.yml/badge.svg?branch=main
[gh-actions-url]: https://github.com/irojkov-ph/CVsim.jl/actions/workflows/CI.yml?query=branch%3Amain

<!-- ![Package Logo](package_logo.png) -->

CVsim.jl is a Julia package for simulating continuous variable quantum systems, with a focus on studying quantum information processing using bosonic codes such as Gottesman-Kitaev-Preskill (GKP), Cat or Binomial states. Unlike traditional approaches that employ the Fock basis representation of states, this package utilizes the first quantization method, specifically sampling the position or momentum quadratures.
 
## Features

- Simulation of continuous variable quantum systems
- Utilizes the first quantization approach for state and operator representation
- Emphasis on bosonic codes (e.g. Cat, GKP and Binomial codes), including some theoretically suggested and/or experimentally demonstrated protocols for:
  - State preparation 
  - Error-correction and stabilization
- Fast calculation of the Wigner quasiprobability function using fast Fourier transfrom (FFT)


## Installation

To use this package, you need to have Julia installed on your system. You can download Julia from the official website: [https://julialang.org/](https://julialang.org/)

CVsim.jl should be available in the Julia General Registry such that its installation can simply be done by
```
   using Pkg; Pkg.add("CVsim")
```
or through the package manager:
1. Open a Julia REPL (Read-Eval-Print Loop) or a Julia IDE.
2. Change the directory to the root of the cloned/downloaded repository.
3. Enter the package manager by pressing `]`.
4. Activate the package by running the following command:
   ```
    add CVsim
   ```
   or alternatively (if you want the latest, potentially unstable, release)
   ```
    add https://github.com/irojkov-ph/CVsim.jl
   ```
7. After installation, exit the package manager by pressing `backspace`.

## Usage

To use the package, follow these steps:

1. Open a Julia REPL or a Julia IDE.
2. Import the package:
   ```julia
    using CVsim
   ```
3. Create and initialize the parameters of a continuous variable quantum system:
   ```julia
      P = CVsim_Parameters(qvec_min=..., qvec_max=..., n_points=..., b_prefered=...)
   ```
   - `qvec_min` (Float): Minimum position.
   - `qvec_max` (Float): Maximum position.
   - `n_points` (Int): Number of points equidistantly distributed between `qvec_min` and `qvec_max`.
   - `encoding` (String): Preferred basis for state representation which is either `"position"` or `"momentum"`.
4. Perform desired operations and calculations using the provided functions and methods.

For detailed documentation, please refer to the [documentation folder](./docs) in this repository.

## Examples

### Simulating a GKP State

```julia
  using CVsim

  # Create a continuous variable quantum system with default parameters
  system = CVsim_Parameters()

  # Initialize a default GKP state (which is a logical |0⟩ state)
  GKP0 = prepare_GKP(P)

  # Perform operations on the state
  # ...

  # Calculate the Wigner quasi probability function of the state
  Wig = wignerify(GKP0)

  # Plot the result (here we use the Plots.jl package for plotting)
  using Plots; pyplot();

  ticks  = [ k*sqrt(pi) for k = div(P.b.xmin,sqrt(pi)):2:div(P.b.xmax,sqrt(pi)) ];
  labels = [ string(round(k)) for k = div(P.b.xmin,sqrt(pi)):2:div(P.b.xmax,sqrt(pi)) ];

  contourf(P.qvec, P.pvec, Wig, legend=false,
           c=cgrad([:blue, :white, :red]),
           clim=(-maximum(Wig),maximum(Wig)),
           levels=200, framestyle = :box, aspect_ratio=1,
           alpha=0.7,grid=true,gridalpha=1.0,
           ylabel="p", ylim = (ticks[1],ticks[end]),
           yticks=((ticks,labels)),
           xlabel="q", xlim = (ticks[1],ticks[end]), 
           xticks=((ticks,labels))
          )
```

## Citing

This package was developed and employed for the first time for studying gates between two finite-energy GKP states (see [Rojkov et al. (2023)](https://doi.org/10.48550/arXiv.2305.05262)). If you use this package in the context of a potential publication we would be pleased if you cite our work as well as the package itself using the following bibtex entries

```bibtex
  @article{two_qubit_rojkov_2023,
    title     = {Two-qubit operations for finite-energy Gottesman-Kitaev-Preskill encodings}, 
    author    = {Ivan Rojkov and Paul Moser R\"{o}ggla and Martin Wagener and Moritz Fontbot\'{e}-Schmidt and Stephan Welte and Jonathan Home and Florentin Reiter},
    journal   = {arXiv:quant-ph/2305.05262},
    year      = {2023},
    month     = may,
  }
```

```bibtex
  @software{cvsim_rojkov_2024,
    author    = {Rojkov, Ivan},
    title     = {CVsim.jl: Continuous variable quantum system simulation package},
    month     = jul,
    year      = 2024,
    publisher = {Zenodo},
    version   = {v0.1.0},
    doi       = {10.5281/zenodo.13127577},
    url       = {https://doi.org/10.5281/zenodo.13127577}
  }
```

## Contributing

Contributions to this package are welcome! If you encounter any issues or have suggestions for improvements, please create a new issue in the [issue tracker](https://github.com/irojkov-ph/CVsim.jl/issues).

If you would like to contribute code, fork the repository, make your changes, and submit a pull request. Please ensure that your code follows the established coding style and is well-documented.

## License

This package is distributed under the MIT License.