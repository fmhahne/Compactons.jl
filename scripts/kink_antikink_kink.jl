using Compactons
using DrWatson
include(srcdir("plots.jl"))

v = 0.5
vL = v
vR = -v

params = KinkAntikinkKink(; vL, vR, tmax=20.0)
data, _ = produce_or_load(datadir("kink_antikink_kink"), params, simulation)
@unpack x, t, η, H = data

fig, axs = plt.subplots(1, 2; figsize=(6, 3))
heatmap!(axs[1], x, t, η; colorbar=true)
axs[1].set_xlim(-10, 10)
axs[1].set_title(raw"$\eta(t,x)$")
axs[1].label_outer()

heatmap!(axs[2], x, t, H; colorbar=true, cmap="magma", norm=mpl.colors.SymLogNorm(1e-5))
axs[2].set_xlim(-10, 10)
axs[2].set_title(raw"$\mathcal{H}(t,x)$")
axs[2].label_outer()

Base.mkpath(plotsdir("kink_antikink_kink"))
fig.savefig(plotsdir("kink_antikink_kink", "example.pdf"))
fig
