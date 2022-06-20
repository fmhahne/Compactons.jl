using DrWatson
using Compactons
include(srcdir("plots.jl"))

V = 0.65
params = KinkKink(V)
data, _ = produce_or_load(datadir("kink_kink"), params, simulation)
@unpack x, t, η = data

fig, ax = plt.subplots()
heatmap!(ax, x, t, η; colorbar=true, norm=mpl.colors.CenteredNorm())
ax.set_xlim(-10V, 10V)
fig.savefig(plotsdir("kink_kink", savename(params, "pdf")))
fig
