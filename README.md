# Compactons.jl

Simulations for interacting compactons in scalar field theories with non-analytical (V-shaped) potentials.

## Prequisites

Install official Julia binaries:

```sh
curl -fsSL https://install.julialang.org | sh
```

Install required Julia packages:

```
julia
julia> ]
(@v1.7) pkg> activate .
(Compactons) pkg> instantiate
```

## Generating data and plots

```sh
julia --project -t auto scripts/<script name>.jl
```
