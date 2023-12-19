using DrWatson
using Compactons
include(srcdir("plots.jl"))

α = 0.0
v₀ = 0.0
v = 0.95

params = SymOscillonScattering(; v, v₀, α)
data, _ = produce_or_load(datadir("oscillon_scattering"), params, simulation)
@unpack x, t, ϕ, H = data

fig, ax = plt.subplots()
vmax = maximum(abs.(ϕ))
norm = mpl.colors.SymLogNorm(1e-4; vmin=-vmax, vmax)
heatmap!(ax, x, t, ϕ; cmap="RdBu", norm, colorbar=true)
ax.set_xlim(-3, 3)
ax.set_ylim(0, 4)

mkpath(plotsdir("oscillon_scattering"))
fig.savefig(plotsdir("oscillon_scattering", "example.pdf"))
fig
