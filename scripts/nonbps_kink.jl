using DrWatson
using Compactons
include(srcdir("plots.jl"))

let ϵs = [-0.15, 0.15, -0.30, 0.30]
    fig, axs = plt.subplots(2, 2; figsize=(6.2, 5.8), sharex=true, sharey=true, tight_layout=false)

    cmap = mpl.cm.get_cmap("magma")
    norm = mpl.colors.SymLogNorm(1e-5, clip=true)
    cb = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

    for (ϵ, ax) ∈ zip(ϵs, axs)
        data, _ = produce_or_load(datadir("nonbps_kink"), NonBPSKink(ϵ), simulation)
        @unpack x, t, η, H = data

        ax.imshow(H'; origin="lower", extent=[x[begin], x[end], t[begin], t[end]], cmap=cmap, norm=norm)
        ax.set_xlim(-5, 5)
        ax.set_title("\$\\epsilon = $ϵ \$")

        cb.set_array(H)
        cb.autoscale()
    end

    axs[1, 1].set_ylabel(raw"$t$")
    axs[2, 1].set_ylabel(raw"$t$")
    axs[2, 1].set_xlabel(raw"$x$")
    axs[2, 2].set_xlabel(raw"$x$")

    fig.tight_layout()
    fig.subplots_adjust(right=0.85)
    cax = fig.add_axes([0.87, 0.15, 0.015, 0.7])
    fig.colorbar(cb, cax=cax)

    fig.savefig(plotsdir("nonbps_kink", "hamiltonian.pdf"))
    fig
end
