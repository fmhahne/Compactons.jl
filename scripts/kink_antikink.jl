using DrWatson
using Compactons
include(srcdir("plots.jl"))

fig, axs = plt.subplots(3, 2; figsize=(6.2, 5), sharex=true, sharey=true)

for (i, V) ∈ enumerate([0.3, 0.6, 0.9])
    data, _ = produce_or_load(datadir("kink_antikink"), KinkAntikink(V), simulation)
    @unpack x, t, η, H = data

    heatmap!(axs[i, 1], x, t, η; colorbar=true, norm=mpl.colors.CenteredNorm())
    axs[i, 1].set_xlim(-10, 10)
    axs[i, 1].set_title(raw"$\eta(t,x)\quad V=" * "$V" * raw"$")
    axs[i, 1].label_outer()

    heatmap!(axs[i, 2], x, t, H; colorbar=true, cmap="magma", norm=mpl.colors.SymLogNorm(1e-5))
    axs[i, 2].set_xlim(-10, 10)
    axs[i, 2].set_title(raw"$\mathcal{H}(t,x)\quad V=" * "$V" * raw"$")
    axs[i, 2].label_outer()
end

fig.savefig(plotsdir("kink_antikink", "examples.pdf"))
fig
