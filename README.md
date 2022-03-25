# Compactons.jl

Simulations for interacting compactons in scalar field theories with non-analitycal (V-shaped) potentials.

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

Install required Python dependencies:

```sh
python -m pip install -r requirements.txt
```

## Generating data and figures

```sh
cd <simulation name> # i.e. kink-oscillon-scattering
julia --project=.. -t auto <simulation name>.jl
python fig.py
```
