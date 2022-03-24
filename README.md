# Compactons.jl

## Install oficial Julia binaries

```sh
curl -fsSL https://install.julialang.org | sh
```

## Install Julia dependecies

```julia
julia
julia> ]
(@v1.7) pkg> activate .
(Compactons.jl) pkg> instantiate
```

## Install Python dependencies

```sh
python -m pip install -r requirements.txt
```

## Generate data from simulation

```sh
cd <simulation name> # i.e. kink-oscillon-scattering
mkdir -p data
julia --project=.. -t auto <simulation name>.jl # i.e. kink-oscillon-scattering.jl
```

## Generate heatmaps from data

```sh
cd <simulation name> # i.e. kink-oscillon-scattering
mkdir -p fig
python fig.py
```
