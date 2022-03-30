# Compactons.jl

Simulations for interacting compactons in scalar field theories with non-analytical (V-shaped) potentials.

## Prequisites

Install oficial Julia binaries:

```sh
curl -fsSL https://install.julialang.org | sh
```

Install required Julia packages:

```julia
julia
julia> ]
(@v1.7) pkg> activate .
(Compactons.jl) pkg> instantiate
```

Install required Python packages:

```sh
python -m pip install -r requirements.txt
```

## Generating data and figures

```sh
julia --project -t auto scripts/<simulation name>.jl
python scripts/<simulation name>_plots.py
```
