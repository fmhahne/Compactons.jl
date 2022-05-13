using DrWatson
using Compactons
include(srcdir("plots.jl"))

xlim = (-5, 5) .+ π / 2

let V = 0.0, v₀ = 0.0, l = 1.0
    fig, axs = plt.subplots(1, 2; figsize=(6.2, 3.1), sharey=true, tight_layout=false)

    norm = mpl.colors.SymLogNorm(1e-5, clip=true)
    cmap = mpl.cm.get_cmap("magma")
    cb = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

    for (α, ax) ∈ zip([0, 0.25], axs)
        x₀ = -(π - x_R(α, V; l, v₀) - x_L(α, V; l, v₀)) / 2
        data, _ = produce_or_load(datadir("kink_oscillon_superposition"), KinkOscillon(; l, V, α, v₀, x₀), simulation)
        @unpack x, t, η, H = data

        ax.imshow(H'; origin="lower", extent=[x[begin], x[end], t[begin], t[end]], cmap=cmap, norm=norm)

        ax.set_xlim(xlim...)
        ax.set_title("\$\\alpha = $α\$")

        ax.set_xlabel(raw"$x$")

        cb.set_array(H)
        cb.autoscale()
    end

    axs[1].set_ylabel(raw"$t$")

    fig.tight_layout()
    fig.subplots_adjust(right=0.85)
    cax = fig.add_axes([0.87, 0.16, 0.015, 0.74])
    fig.colorbar(cb, cax=cax)

    fig.savefig(plotsdir("kink_oscillon_superposition", "non_zero_p.pdf"))
    fig
end

let V = 0.0, v₀ = 0.0, l = 1.0
    fig, axs = plt.subplots(1, 2; figsize=(6.2, 3.1), sharey=true, tight_layout=false)

    norm = mpl.colors.CenteredNorm()
    cmap = mpl.cm.get_cmap("RdBu")
    cb = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

    for (α, ax) ∈ zip([0.25, 0.75], axs)
        x₀ = -(π - x_R(α, V; l, v₀) - x_L(α, V; l, v₀)) / 2
        data, _ = produce_or_load(datadir("kink_oscillon_superposition"), KinkOscillon(; l, V, α, v₀, x₀), simulation)
        @unpack x, t, η, H = data

        χ = η - kink.(t', x)
        ax.imshow(χ'; origin="lower", extent=[x[begin], x[end], t[begin], t[end]], cmap=cmap, norm=norm)

        ax.set_xlim(xlim...)

        ax.set_title("\$\\alpha = $α\$")
        ax.set_xlabel(raw"$x$")

        ax.axvline(0; linewidth=0.5, color="black", linestyle="dashed")
        ax.axvline(float(π); linewidth=0.5, color="black", linestyle="dashed")

        cb.set_array(χ)
        cb.autoscale()
    end

    axs[1].set_ylabel(raw"$t$")

    fig.tight_layout()
    fig.subplots_adjust(right=0.85)
    cax = fig.add_axes([0.87, 0.16, 0.015, 0.74])
    fig.colorbar(cb, cax=cax)

    fig.savefig(plotsdir("kink_oscillon_superposition", "zero_p.pdf"))
    fig
end

let l = 1.0, V = 0.75, v₀ = 0.0
    fig, axs = plt.subplots(2, 2; figsize=(6.2, 5.8), sharex=true, sharey=true, tight_layout=true)

    for (i, α) ∈ enumerate([0.0, 0.25])
        x₀ = -(π - x_R(α, V; l, v₀) - x_L(α, V; l, v₀)) / 2
        data, _ = produce_or_load(datadir("kink_oscillon_superposition"), KinkOscillon(; l, V, α, v₀, x₀), simulation)
        @unpack x, t, η, H = data

        χ = η - kink.(t', x)

        axs[i, 1].set_title("\$\\chi \\, (\\alpha = $α)\$")
        heatmap!(axs[i, 1], x, t, χ; cmap="RdBu", norm=mpl.colors.CenteredNorm())

        axs[i, 1].set_xlim(-2, 8)
        axs[i, 1].set_xlabel(nothing)
        axs[i, 1].set_ylabel(nothing)
        axs[i, 1].axvline(0; linewidth=0.5, color="black", linestyle="dashed")
        axs[i, 1].axvline(float(π); linewidth=0.5, color="black", linestyle="dashed")

        axs[i, 2].set_title("\$\\mathcal{H} \\, (\\alpha = $α)\$")
        heatmap!(axs[i, 2], x, t, H; cmap="magma", norm=mpl.colors.SymLogNorm(1e-5))

        axs[i, 2].set_xlim(-2, 8)
        axs[i, 2].set_xlabel(nothing)
        axs[i, 2].set_ylabel(nothing)
        axs[i, 2].axvline(0; linewidth=0.5, color="lime", linestyle="dashed")
        axs[i, 2].axvline(float(π); linewidth=0.5, color="lime", linestyle="dashed")
    end

    axs[1, 1].set_ylabel(raw"$t$")
    axs[2, 1].set_ylabel(raw"$t$")
    axs[2, 1].set_xlabel(raw"$x$")
    axs[2, 2].set_xlabel(raw"$x$")

    fig.savefig(plotsdir("kink_oscillon_superposition", "centered-alpha.pdf"))
    fig
end

let l = 1.0, V = 0.75, v₀ = 0.0
    fig, axs = plt.subplots(2, 2; figsize=(6.2, 5.8), sharex=true, sharey=true, tight_layout=true)

    for (i, α) ∈ enumerate([0.0, 0.25])
        x₀ = x_L(α, V; l, v₀)
        data, _ = produce_or_load(datadir("kink_oscillon_superposition"), KinkOscillon(; l, V, α, v₀, x₀), simulation)
        @unpack x, t, η, H = data

        χ = η - kink.(t', x)

        axs[i, 1].set_title("\$\\chi \\, (\\alpha = $α)\$")
        heatmap!(axs[i, 1], x, t, χ; cmap="RdBu", norm=mpl.colors.CenteredNorm())

        axs[i, 1].set_xlim(-2, 8)
        axs[i, 1].set_xlabel(nothing)
        axs[i, 1].set_ylabel(nothing)
        axs[i, 1].axvline(0; linewidth=0.5, color="black", linestyle="dashed")
        axs[i, 1].axvline(float(π); linewidth=0.5, color="black", linestyle="dashed")

        axs[i, 2].set_title("\$\\mathcal{H} \\, (\\alpha = $α)\$")
        heatmap!(axs[i, 2], x, t, H; cmap="magma", norm=mpl.colors.SymLogNorm(1e-5))

        axs[i, 2].set_xlim(-2, 8)
        axs[i, 2].set_xlabel(nothing)
        axs[i, 2].set_ylabel(nothing)
        axs[i, 2].axvline(0; linewidth=0.5, color="lime", linestyle="dashed")
        axs[i, 2].axvline(float(π); linewidth=0.5, color="lime", linestyle="dashed")
    end

    axs[1, 1].set_ylabel(raw"$t$")
    axs[2, 1].set_ylabel(raw"$t$")
    axs[2, 1].set_xlabel(raw"$x$")
    axs[2, 2].set_xlabel(raw"$x$")

    fig.savefig(plotsdir("kink_oscillon_superposition", "left-alpha.pdf"))
    fig
end
