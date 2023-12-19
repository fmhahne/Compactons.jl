using Compactons
using DrWatson
include(srcdir("plots.jl"))

params = TriangularDelta(; A=0.2, ϵ=0.1)
data, _ = produce_or_load(datadir("triangular_delta"), params, simulation)
@unpack x, t, ϕ, H = data

fig, ax = plt.subplots()
vmax = maximum(abs.(ϕ))
norm = mpl.colors.SymLogNorm(1e-2; vmin=-vmax, vmax)
heatmap!(ax, x, t, ϕ; colorbar=true, cmap="RdBu", norm)
ax.set_xlim(-6, 6)
Base.mkpath(plotsdir("triangular_delta"))
fig.savefig(plotsdir("triangular_delta", "example.pdf"))
fig
