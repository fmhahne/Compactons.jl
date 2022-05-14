using DrWatson
using Compactons
include(srcdir("plots.jl"))

ls = 0.5:0.5:2.0
Vs = 0.0:0.25:0.75
αs = 0.0:0.25:0.75
v₀s = 0.0:0.25:0.75

xlim = (-2, 8)

let V = 0.75, α = 0.0, v₀ = 0.0
    fig, axs = plt.subplots(2, 2; figsize=(6.2, 5.5), sharex=true, sharey=true)

    for (l, ax) ∈ zip(ls, Iterators.flatten(eachrow(axs)))
        data, _ = produce_or_load(datadir("kink_oscillon_scattering"), KinkOscillon(; l, V, α, v₀), simulation)
        @unpack x, t, η, H, E₁, E₂, E₃ = data

        χ = η - kink.(t', x)
        heatmap!(ax, x, t, χ; colorbar=true, cmap="RdBu", norm=mpl.colors.CenteredNorm())
        show_kink_borders!(ax)

        ax.set_xlim(xlim...)
        ax.set_title("\$l = $l\$")
        ax.label_outer()
    end

    fig.savefig(plotsdir("kink_oscillon_scattering", "l.pdf"))
    fig
end

let l = 0.5, α = 0.0, v₀ = 0.0
    fig, axs = plt.subplots(2, 2; figsize=(6.2, 5.8), sharex=true, sharey=true, tight_layout=false)

    norm = mpl.colors.CenteredNorm()
    cmap = mpl.cm.get_cmap("RdBu")
    cb = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

    for (V, ax) ∈ zip(Vs, Iterators.flatten(eachrow(axs)))
        data, _ = produce_or_load(datadir("kink_oscillon_scattering"), KinkOscillon(; l, V, α, v₀), simulation)
        @unpack x, t, η, H, E₁, E₂, E₃ = data

        χ = η - kink.(t', x)
        heatmap!(ax, x, t, χ; cmap="RdBu", norm=mpl.colors.CenteredNorm())
        show_kink_borders!(ax)

        ax.set_xlim(xlim...)
        ax.set_title("\$V = $V\$")
        ax.label_outer()

        cb.set_array(χ)
        cb.autoscale()
    end

    shared_colorbar!(fig, cb)

    fig.savefig(plotsdir("kink_oscillon_scattering", "V.pdf"))
    fig
end

let l = 1.0, V = 0.6, v₀ = 0.0
    fig, axs = plt.subplots(2, 2; figsize=(6.2, 5.8), sharex=true, sharey=true, tight_layout=false)

    norm = mpl.colors.CenteredNorm()
    cmap = mpl.cm.get_cmap("RdBu")
    cb = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

    for (α, ax) ∈ zip(αs, Iterators.flatten(eachrow(axs)))
        data, _ = produce_or_load(datadir("kink_oscillon_scattering"), KinkOscillon(; l, V, α, v₀), simulation)
        @unpack x, t, η, H, E₁, E₂, E₃ = data

        χ = η - kink.(t', x)
        heatmap!(ax, x, t, χ; cmap="RdBu", norm=mpl.colors.CenteredNorm())
        show_kink_borders!(ax)

        ax.set_xlim(xlim...)
        ax.set_title("\$\\alpha = $α\$")
        ax.label_outer()

        cb.set_array(χ)
        cb.autoscale()
    end

    shared_colorbar!(fig, cb)

    fig.savefig(plotsdir("kink_oscillon_scattering", "alpha.pdf"))
    fig
end

let l = 0.75, V = 0.8, α = 0.0
    fig, axs = plt.subplots(2, 2; figsize=(6.2, 5.8), sharex=true, sharey=true, tight_layout=false)

    norm = mpl.colors.CenteredNorm()
    cmap = mpl.cm.get_cmap("RdBu")
    cb = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

    for (v₀, ax) ∈ zip(v₀s, Iterators.flatten(eachrow(axs)))
        data, _ = produce_or_load(datadir("kink_oscillon_scattering"), KinkOscillon(; l, V, α, v₀), simulation)
        @unpack x, t, η, H, E₁, E₂, E₃ = data

        χ = η - kink.(t', x)
        heatmap!(ax, x, t, χ; cmap="RdBu", norm=mpl.colors.CenteredNorm())
        show_kink_borders!(ax)

        ax.set_xlim(xlim...)
        ax.set_title("\$v_0 = $v₀\$")
        ax.label_outer()

        cb.set_array(χ)
        cb.autoscale()
    end

    shared_colorbar!(fig, cb)

    fig.savefig(plotsdir("kink_oscillon_scattering", "v0.pdf"))
    fig
end

let l = 1.0, Vs = [0.6, 0.75], αs = [0.0, 0.50], v₀ = 0.0
    fig, axs = plt.subplots(2, 2; figsize=(6.2, 6.2 * 2 / (1 + √5)), sharex=true, sharey=true, tight_layout=false)

    for ((V, α), ax) ∈ zip(Iterators.product(Vs, αs), Iterators.flatten(eachrow(axs)))
        data, _ = produce_or_load(datadir("kink_oscillon_scattering"), KinkOscillon(; l, V, α, v₀), simulation)
        @unpack x, t, η, H, E₁, E₂, E₃ = data

        ax.plot(t, E₁ / E₁[begin]; label=raw"$E_1(t) / E_\mathrm{osc}$")
        ax.plot(t, E₂ / E₁[begin]; label=raw"$E_2(t) / E_\mathrm{osc}$")
        ax.plot(t, E₃ / E₁[begin]; label=raw"$E_3(t) / E_\mathrm{osc}$")

        ax.set_title("\$V = $V\$, \$\\alpha = $α\$")
        ax.label_outer()
    end

    fig.tight_layout()
    handles, labels = axs[1, 1].get_legend_handles_labels()
    fig.subplots_adjust(bottom=0.15)
    fig.legend(handles, labels; loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0))

    fig.savefig(plotsdir("kink_oscillon_scattering", "energies-vs-time.pdf"))
    fig
end

let l = 0.5, V = 0.8, α = 0.0, v₀ = 0.0
    fig, ax = plt.subplots(sharex=true, sharey=true)

    data, _ = produce_or_load(datadir("kink_oscillon_scattering"), KinkOscillon(; l, V, α, v₀), simulation)
    @unpack x, t, η, H, E₁, E₂, E₃ = data

    χ = η - kink.(t', x)

    ax.plot(x, χ[:, begin]; label=raw"$ t = 0 $")

    n = findfirst(t .== 0.15)
    ax.plot(x, χ[:, n]; label="\$ t = $(t[n]) \$")

    n = findfirst(t .== 0.30)
    ax.plot(x, χ[:, n]; label="\$ t = $(t[n]) \$")

    show_kink_borders!(ax)

    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel(raw"$x$")

    ax.legend()

    fig.savefig(plotsdir("kink_oscillon_scattering", "profile.pdf"))
    fig
end
