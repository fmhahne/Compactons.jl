# Compactons.jl

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

Simulations for interacting compactons in scalar field theories with non-analytical (V-shaped) potentials.

The code on this repository was developed for use in the following publications:

- Hahne, F.M., Klimas, P. [Compact kink and its interaction with compact oscillons](https://doi.org/10.1007/JHEP09(2022)100). *J. High Energ. Phys.* 2022, 100 (2022)
- Hahne, F.M., Klimas, P. [Scattering of compact kinks](https://doi.org/10.1007/JHEP01(2024)067). *J. High Energ. Phys.* 2024, 67 (2024)

## Prerequisites

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
