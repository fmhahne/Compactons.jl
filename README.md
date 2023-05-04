# KinkOscillonInteraction.jl

Julia module and scripts used to produce the simulations for the paper "[Compact kink and its interaction with compact oscillons](https://arxiv.org/abs/2207.07064)".

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
(KinkOscillonInteraction) pkg> instantiate
```

## Generating data and plots

```sh
julia --project -t auto scripts/<script name>.jl
```
