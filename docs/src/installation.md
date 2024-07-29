# Installation

To use this package, you need to have Julia installed on your system. You can download Julia from the official website: [https://julialang.org/](https://julialang.org/)

CVsim.jl should be available in the Julia General Registry such that its installation can simply be done by
```julia
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

5. After installation, exit the package manager by pressing `backspace`.
