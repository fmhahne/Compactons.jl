using DrWatson
using DifferentialEquations
using Compactons
include(srcdir("plots.jl"))

let l = 1.0, V = 0.0, α = 0.25, v₀ = 0.0, model = toy
    x₀ = -(2√2 - x_R(α, V; l, v₀) - x_L(α, V; l, v₀)) / 2
    parameters = KinkOscillon(; l, V, α, v₀, x₀, model)
    sim(p) = simulation(p; dx=5e-4, sampling=20)
    data, _ = produce_or_load(datadir("toy_model"), parameters, sim)
    @unpack x, t, η, H = data
    χ = η - toykink.(t', x)

    fig, axs = plt.subplots(1, 2; figsize=(6.2, 2.8))
    heatmap!(axs[1, 1], x, t, χ; norm=mpl.colors.CenteredNorm(), cmap="RdBu", colorbar=true)
    axs[1].set_xlim(√2 - 5, √2 + 5)
    axs[1].axvline(0; linewidth=0.5, color="black", linestyle="dashed")
    axs[1].axvline(2√2; linewidth=0.5, color="black", linestyle="dashed")
    axs[1].set_title(raw"$\chi(t, x)$")

    n = 76
    axs[2].plot(x, χ[:, n]; color="black")
    axs[2].set_xlim(0, 2√2)
    axs[2].set_xlabel("x")
    axs[2].set_title("\$ \\chi($(t[n]), x) \$")

    fig.savefig(plotsdir("toy_model", "centered.pdf"))
    fig
end

let l = 0.2, V = 0.0, α = 0.25, v₀ = 0.0, model = toy
    x₀ = -(√2 - x_R(α, V; l, v₀) - x_L(α, V; l, v₀)) / 2
    parameters = KinkOscillon(; l, V, α, v₀, x₀, model)
    sim(p) = simulation(p; dx=5e-4, sampling=20)
    data, _ = produce_or_load(datadir("toy_model"), parameters, sim)
    @unpack x, t, η, H = data
    χ = η - toykink.(t', x)

    fig, axs = plt.subplots(1, 2; figsize=(6.2, 2.8))
    heatmap!(axs[1, 1], x, t, χ; norm=mpl.colors.CenteredNorm(), cmap="RdBu", colorbar=true)
    axs[1].set_xlim(√2 - 5, √2 + 5)
    axs[1].axvline(0; linewidth=0.5, color="black", linestyle="dashed")
    axs[1].axvline(2√2; linewidth=0.5, color="black", linestyle="dashed")
    axs[1].set_title(raw"$\chi(t, x)$")

    n = 21
    axs[2].plot(x, χ[:, n]; color="black")
    axs[2].set_xlim(0, √2)
    axs[2].set_xlabel("x")
    axs[2].set_title("\$ \\chi($(t[n]), x) \$")

    fig.savefig(plotsdir("toy_model", "quarter-rest.pdf"))
    fig
end

let l = 0.2, V = 0.75, α = 0.0, v₀ = 0.0, model = toy
    x₀ = -(√2 - x_R(α, V; l, v₀) - x_L(α, V; l, v₀)) / 2
    parameters = KinkOscillon(; l, V, α, v₀, x₀, model)
    sim(p) = simulation(p; dx=5e-4, sampling=20)
    data, _ = produce_or_load(datadir("toy_model"), parameters, sim)
    @unpack x, t, η, H = data
    χ = η - toykink.(t', x)

    fig, axs = plt.subplots(1, 2; figsize=(6.2, 2.8))
    heatmap!(axs[1, 1], x, t, χ; norm=mpl.colors.CenteredNorm(), cmap="RdBu", colorbar=true)
    axs[1].set_xlim(√2 - 5, √2 + 5)
    axs[1].axvline(0; linewidth=0.5, color="black", linestyle="dashed")
    axs[1].axvline(2√2; linewidth=0.5, color="black", linestyle="dashed")
    axs[1].set_title(raw"$\chi(t, x)$")

    n = 21
    axs[2].plot(x, χ[:, n]; color="black")
    axs[2].set_xlim(0, √2)
    axs[2].set_xlabel("x")
    axs[2].set_title("\$ \\chi($(t[n]), x) \$")

    fig.savefig(plotsdir("toy_model", "quarter-moving.pdf"))
    fig
end
