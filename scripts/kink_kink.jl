using DrWatson
using Compactons
include(srcdir("plots.jl"))

V = 0.65
data, _ = produce_or_load(datadir("kink_kink"), KinkKink(V), simulation)
@unpack x, t, η, H = data

fig, axs = plt.subplots(1, 2; sharey=true)
heatmap!(axs[1], x, t, η; colorbar=true, norm=mpl.colors.CenteredNorm())
heatmap!(axs[2], x, t, H; colorbar=true, cmap="magma", norm=mpl.colors.SymLogNorm(1e-5))
axs[2].set_ylabel(nothing)
fig
