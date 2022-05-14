using DrWatson
using Compactons
include(srcdir("plots.jl"))

let ϵs = [-0.15, 0.15, -0.30, 0.30]
    fig, axs = plt.subplots(2, 2; figsize=(6.2, 5.8), sharex=true, sharey=true, tight_layout=false)

    cmap = mpl.cm.get_cmap("magma")
    norm = mpl.colors.SymLogNorm(1e-5, clip=true)
    cb = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

    for (ϵ, ax) ∈ zip(ϵs, Iterators.flatten(eachrow(axs)))
        data, _ = produce_or_load(datadir("nonbps_kink"), NonBPSKink(ϵ), simulation)
        @unpack x, t, η, H = data

        heatmap!(ax, x, t, H; cmap=cmap, norm=norm)

        ax.set_xlim(-5, 5)
        ax.set_title("\$\\epsilon = $ϵ \$")
        ax.label_outer()

        cb.set_array(H)
        cb.autoscale()
    end

    shared_colorbar!(fig, cb)

    fig.savefig(plotsdir("nonbps_kink", "hamiltonian.pdf"))
    fig
end
