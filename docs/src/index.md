# CVsim.jl: Continuous variable quantum system simulation package 

`CVsim.jl` is a Julia package for simulating continuous variable (CV) quantum systems, with a focus on studying quantum information processing using bosonic codes such as Gottesman-Kitaev-Preskill (GKP), Cat or Binomial states. Unlike traditional approaches that employ the Fock basis representation of states, this package utilizes the first quantization method, specifically sampling the position or momentum quadratures.

## Features

The package allows you to setup the settings of a CV as well as of an auxiliary discrete variable (qubit or qudit) quantum systems, prepare and stabilize some CV codewords using different methods proposed and/or demonstrated in the last decades, and calculate the Wigner quasiprobability distribution of these states. All these elements are presented in detail in the corresponding sections of the documentation.

## Current limitations

So far mainly the preparation and stabilization of GKP states have been implemented, but we believe that representing other continuous variable states in the position or momentum bases can help in some cases alleviate the computational complexity associated with using the Fock basis. 

## Potential improvements

We believe that using position and momentum representations of quantum states and operators can be provide an efficient way to simulate CV system compared of more standard strategies (especially because of the Fourier transform relation between the two bases). We are currently not planning any publication regarding the performance of CVsim.jl. However, if you like the package, please contact us, we will be happy to collaborate!

If you are working with CV encodings, we will be happy of your contributions to the package. Here are some potential improvements that could be envisaged:
- State preparation and stabilization schemes for Cat and Binomial qubits, e.g. [Mirrahimi _et al._ (2014)](https://iopscience.iop.org/article/10.1088/1367-2630/16/4/045014) and [Lihm _et al._ (2018)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.98.012317),
- Additional state preparation and stabilization schemes for GKP qubits, especially those developed for photonic platforms such as [Glancy and Knill (2006)](https://doi.org/10.1103/PhysRevA.73.012325) and [Walshe _et al._ (2020)](https://doi.org/10.1103/PhysRevA.102.062411),
- Logical operations and particularly the bias preserving gates for Cat qubits, see [Puri _et al._ (2020)](https://doi.org/10.1126/sciadv.aay5901),
- Experimentally realistic readout protocols for CV state, e.g. [Rosenblum, Reinhold _et al._ (2018)](https://doi.org/10.1126/science.aat3996) for Cat qubits and [Hastrup and Andersen (2021)](https://doi.org/10.1088/2058-9565/ac070d) for GKP codes,
- Joint (i.e. two-mode) Wigner function quasiprobability function calculation, see [Wang _et al._ (2016)](https://doi.org/10.1126/science.aaf2941),
- Characteristic function calculation,
- Dominant noise processes for continuous and discrete variable systems,
- Proper compatibility with `Float32`,
- Unit tests for the package.

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
