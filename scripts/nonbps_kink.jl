using DelimitedFiles
using DrWatson
using Compactons
include(srcdir("plots.jl"))

let ϵs = [-0.15, 0.15, -0.30, 0.30]
    fig, axs = plt.subplots(
        2, 2; figsize=(6.2, 5.8), sharex=true, sharey=true, tight_layout=false
    )

    cmap = mpl.cm.get_cmap("magma")
    norm = mpl.colors.SymLogNorm(1e-5; clip=true)
    cb = mpl.cm.ScalarMappable(; norm=norm, cmap=cmap)

    for (ϵ, ax) in zip(ϵs, Iterators.flatten(eachrow(axs)))
        data, _ = produce_or_load(datadir("nonbps_kink"), NonBPSKink(ϵ), simulation)
        @unpack x, t, η, H = data

        b = readdlm(datadir("perturbed_kink", "borders", "$ϵ.csv"), ',', Float64, '\n')

        heatmap!(ax, x, t, H; cmap=cmap, norm=norm)
        ax.plot(π ./ (2 * b[:, 2]), border[:, 1]; color="lime")
        ax.plot(-π ./ (2 * b[:, 2]), border[:, 1]; color="lime")

        ax.set_xlim(-5, 5)
        ax.set_title("\$\\epsilon = $ϵ \$")
        ax.label_outer()

        cb.set_array(H)
        cb.autoscale()
    end

    shared_colorbar!(fig, cb)

    fig.savefig(plotsdir("nonbps_kink", "hamiltonian-borders.pdf"))
    fig
end

let ϵ = -0.2
    fig, axs = plt.subplots(1, 2; figsize=(6.2, 2.9))
    data, _ = produce_or_load(datadir("nonbps_kink"), NonBPSKink(ϵ)) do parameters
        return simulation(parameters; sampling=5)
    end
    @unpack x, t, η, H = data

    heatmap!(axs[1], x, t, H; colorbar=true, cmap="magma", norm=mpl.colors.SymLogNorm(1e-5))
    axs[1].set_xlim(-5, 5)
    axs[1].add_patch(
        mpl.patches.Rectangle(
            (-3.2, 7.5), 1.5, 1.5; facecolor="none", edgecolor="lime", linewidth=0.5
        ),
    )
    axs[1].set_title(raw"$\mathcal{H}(t, x)$")

    nx₁ = findfirst(x .≥ -3.2)
    nx₂ = findfirst(x .≥ -1.7)
    nt₁ = findfirst(t .≥ 7.5)
    nt₂ = findfirst(t .≥ 9)
    heatmap!(
        axs[2],
        x[nx₁:nx₂],
        t[nt₁:nt₂],
        η[nx₁:nx₂, nt₁:nt₂];
        colorbar=true,
        cmap="RdBu",
        norm=mpl.colors.CenteredNorm(),
    )
    axs[2].set_title(raw"$\eta(t, x)$")
    fig.savefig(plotsdir("nonbps_kink", "zoom.pdf"))
    fig
end
