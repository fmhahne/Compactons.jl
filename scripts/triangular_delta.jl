using Compactons
using DrWatson
include(srcdir("plots.jl"))

fig, axs = plt.subplots(3, 2; figsize=(6, 5), sharex="col", sharey="row")

for (i, k) in enumerate([0.0, 1.0, 10.0])
    params = TriangularDelta(; ϵ=1e-2, k)
    data, _ = produce_or_load(datadir("triangular_delta"), params, simulation)
    @unpack x, t, ϕ, H = data

    heatmap!(axs[i, 1], x, t, ϕ; colorbar=true, cmap="RdBu", norm=mpl.colors.CenteredNorm())
    axs[i, 1].set_title("\$\\phi(t, x)\$, \$k=$k\$")
    axs[i, 1].label_outer()

    heatmap!(
        axs[i, 2], x, t, H; colorbar=true, cmap="magma", norm=mpl.colors.SymLogNorm(1e-5)
    )
    axs[i, 2].set_title("\$\\mathcal{H}(t, x)\$, \$k=$k\$")
    axs[i, 2].label_outer()
end

Base.mkpath(plotsdir("triangular_delta"))
fig.savefig(plotsdir("triangular_delta", "example.pdf"))
fig
