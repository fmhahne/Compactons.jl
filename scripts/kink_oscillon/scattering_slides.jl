# Some example scattering plots for presentation at the SBF Spring Meeting 2023

using DrWatson
using Compactons
include(srcdir("plots.jl"))

xlim = (-2, 8)

let l = 1.0, V = 0.25, α = 0.5, v₀ = 0.2
    fig, axs = plt.subplots(
        1, 2; figsize=(6.5, 3), sharex=true, sharey=true, layout="tight"
    )

    params = KinkOscillon(; l, V, α, v₀)
    data, _ = produce_or_load(datadir("kink_oscillon_scattering"), params, simulation)
    @unpack x, t, η, H, E₁, E₂, E₃ = data

    heatmap!(axs[1], x, t, H; colorbar=true, cmap="magma", norm=mpl.colors.SymLogNorm(1e-5))
    # show_kink_borders!(ax)

    axs[1].set_xlim(xlim...)
    axs[1].label_outer()

    χ = η - kink.(t', x)
    heatmap!(axs[2], x, t, χ; colorbar=true, cmap="RdBu", norm=mpl.colors.CenteredNorm())
    show_kink_borders!(axs[2])

    axs[2].set_xlim(xlim...)
    axs[2].label_outer()

    Base.mkpath(plotsdir("kink_oscillon_slides"))
    fig.savefig(plotsdir("kink_oscillon_slides", savename(params, "pdf")))
    fig
end

let l = 1.0, V = 0.75, α = 0.0, v₀ = 0.0
    fig, axs = plt.subplots(
        1, 2; figsize=(6.5, 3), sharex=true, sharey=true, layout="tight"
    )

    params = KinkOscillon(; l, V, α, v₀)
    data, _ = produce_or_load(datadir("kink_oscillon_scattering"), params, simulation)
    @unpack x, t, η, H, E₁, E₂, E₃ = data

    heatmap!(axs[1], x, t, H; colorbar=true, cmap="magma", norm=mpl.colors.SymLogNorm(1e-5))
    # show_kink_borders!(ax)

    axs[1].set_xlim(xlim...)
    axs[1].label_outer()

    χ = η - kink.(t', x)
    heatmap!(axs[2], x, t, χ; colorbar=true, cmap="RdBu", norm=mpl.colors.CenteredNorm())
    show_kink_borders!(axs[2])

    axs[2].set_xlim(xlim...)
    axs[2].label_outer()

    Base.mkpath(plotsdir("kink_oscillon_slides"))
    fig.savefig(plotsdir("kink_oscillon_slides", savename(params, "pdf")))
    fig
end
