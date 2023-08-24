using DrWatson
using Compactons
include(srcdir("plots.jl"))

xlim = (-5, 5) .+ π / 2

let V = 0.0, v₀ = 0.0, l = 1.0
    fig, axs = plt.subplots(1, 2; figsize=(6.2, 3.1), sharey=true, tight_layout=false)

    norm = mpl.colors.SymLogNorm(1e-5; clip=true)
    cmap = mpl.cm.get_cmap("magma")
    cb = mpl.cm.ScalarMappable(; norm=norm, cmap=cmap)

    for (α, ax) in zip([0, 0.25], axs)
        x₀ = -(π - x_R(α, V; l, v₀) - x_L(α, V; l, v₀)) / 2
        data, _ = produce_or_load(
            datadir("kink_oscillon_superposition"),
            KinkOscillon(; l, V, α, v₀, x₀),
            simulation,
        )
        @unpack x, t, η, H = data

        heatmap!(ax, x, t, H; cmap=cmap, norm=norm)

        ax.set_xlim(xlim...)
        ax.set_title("\$\\alpha = $α\$")
        ax.label_outer()

        cb.set_array(H)
        cb.autoscale()
    end

    shared_colorbar!(fig, cb)

    fig.savefig(plotsdir("kink_oscillon_superposition", "non_zero_p.pdf"))
    fig
end

let V = 0.0, v₀ = 0.0, l = 1.0
    fig, axs = plt.subplots(1, 2; figsize=(6.2, 3.1), sharey=true, tight_layout=false)

    norm = mpl.colors.CenteredNorm()
    cmap = mpl.cm.get_cmap("RdBu")
    cb = mpl.cm.ScalarMappable(; norm=norm, cmap=cmap)

    for (α, ax) in zip([0.25, 0.75], axs)
        x₀ = -(π - x_R(α, V; l, v₀) - x_L(α, V; l, v₀)) / 2
        data, _ = produce_or_load(
            datadir("kink_oscillon_superposition"),
            KinkOscillon(; l, V, α, v₀, x₀),
            simulation,
        )
        @unpack x, t, η, H = data

        χ = η - kink.(t', x)
        heatmap!(ax, x, t, χ; cmap=cmap, norm=norm)

        show_kink_borders!(ax)
        ax.set_xlim(xlim...)
        ax.set_title("\$\\alpha = $α\$")
        ax.label_outer()

        cb.set_array(χ)
        cb.autoscale()
    end

    shared_colorbar!(fig, cb)

    fig.savefig(plotsdir("kink_oscillon_superposition", "zero_p.pdf"))
    fig
end

let l = 1.0, V = 0.75, v₀ = 0.0
    fig, axs = plt.subplots(
        2, 2; figsize=(6.2, 5.8), sharex=true, sharey=true, tight_layout=true
    )

    for (i, α) in enumerate([0.0, 0.25])
        x₀ = -(π - x_R(α, V; l, v₀) - x_L(α, V; l, v₀)) / 2
        data, _ = produce_or_load(
            datadir("kink_oscillon_superposition"),
            KinkOscillon(; l, V, α, v₀, x₀),
            simulation,
        )
        @unpack x, t, η, H = data

        χ = η - kink.(t', x)

        axs[i, 1].set_title("\$\\chi \\, (\\alpha = $α)\$")
        heatmap!(
            axs[i, 1], x, t, χ; colorbar=true, cmap="RdBu", norm=mpl.colors.CenteredNorm()
        )

        show_kink_borders!(axs[i, 1]; color="black")
        axs[i, 1].set_xlim(-2, 8)
        axs[i, 1].label_outer()

        axs[i, 2].set_title("\$\\mathcal{H} \\, (\\alpha = $α)\$")
        heatmap!(
            axs[i, 2],
            x,
            t,
            H;
            colorbar=true,
            cmap="magma",
            norm=mpl.colors.SymLogNorm(1e-5),
        )

        show_kink_borders!(axs[i, 2]; color="lime")
        axs[i, 2].set_xlim(-2, 8)
        axs[i, 2].label_outer()
    end

    fig.savefig(plotsdir("kink_oscillon_superposition", "centered-alpha.pdf"))
    fig
end

let l = 1.0, V = 0.75, v₀ = 0.0
    fig, axs = plt.subplots(
        2, 2; figsize=(6.2, 5.8), sharex=true, sharey=true, tight_layout=true
    )

    for (i, α) in enumerate([0.0, 0.25])
        x₀ = x_L(α, V; l, v₀)
        data, _ = produce_or_load(
            datadir("kink_oscillon_superposition"),
            KinkOscillon(; l, V, α, v₀, x₀),
            simulation,
        )
        @unpack x, t, η, H = data

        χ = η - kink.(t', x)

        axs[i, 1].set_title("\$\\chi \\, (\\alpha = $α)\$")
        heatmap!(
            axs[i, 1], x, t, χ; colorbar=true, cmap="RdBu", norm=mpl.colors.CenteredNorm()
        )

        show_kink_borders!(axs[i, 1]; color="black")
        axs[i, 1].set_xlim(-2, 8)
        axs[i, 1].label_outer()

        axs[i, 2].set_title("\$\\mathcal{H} \\, (\\alpha = $α)\$")
        heatmap!(
            axs[i, 2],
            x,
            t,
            H;
            colorbar=true,
            cmap="magma",
            norm=mpl.colors.SymLogNorm(1e-5),
        )

        show_kink_borders!(axs[i, 2]; color="lime")
        axs[i, 2].set_xlim(-2, 8)
        axs[i, 2].label_outer()
    end

    fig.savefig(plotsdir("kink_oscillon_superposition", "left-alpha.pdf"))
    fig
end
