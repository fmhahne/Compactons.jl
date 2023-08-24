using DrWatson
using Compactons
include(srcdir("plots.jl"))

V = 0.65
params = KinkKink(; V)
data, _ = produce_or_load(datadir("kink_kink"), params, simulation)
@unpack x, t, η, H = data

fig, axs = plt.subplots(1, 2; figsize=(6, 2.2), sharey=true)

heatmap!(axs[1], x, t, η; colorbar=true, norm=mpl.colors.CenteredNorm())
axs[1].set_xlim(-7, 7)
axs[1].set_title(raw"$\eta(t,x)$")
axs[1].label_outer()

heatmap!(axs[2], x, t, H; colorbar=true, cmap="magma", norm=mpl.colors.SymLogNorm(1e-5))
axs[2].set_xlim(-7, 7)
axs[2].set_title(raw"$\mathcal{H}(t,x)$")
axs[2].label_outer()

fig.savefig(plotsdir("kink_kink", "example.pdf"))
fig
