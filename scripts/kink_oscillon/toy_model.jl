using DrWatson
using DifferentialEquations
using Compactons
include(srcdir("plots.jl"))

let l = 1.0, V = 0.0, α = 0.25, v₀ = 0.0, model = toy
    x₀ = -(2√2 - x_R(α, V; l, v₀) - x_L(α, V; l, v₀)) / 2
    parameters = KinkOscillon(; l, V, α, v₀, x₀, model)
    data, _ = produce_or_load(datadir("toy_model"), parameters) do parameters
        return simulation(parameters; dx=5e-4, sampling=20)
    end
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
    data, _ = produce_or_load(datadir("toy_model"), parameters) do parameters
        return simulation(parameters; dx=5e-4, sampling=20)
    end
    @unpack x, t, η, H = data
    χ = η - toykink.(t', x)

    fig, axs = plt.subplots(1, 2; figsize=(6.2, 2.8))
    heatmap!(axs[1, 1], x, t, χ; norm=mpl.colors.CenteredNorm(), cmap="RdBu", colorbar=true)
    axs[1].set_xlim(√2 - 5, √2 + 5)
    axs[1].axvline(0; linewidth=0.5, color="black", linestyle="dashed")
    axs[1].axvline(2√2; linewidth=0.5, color="black", linestyle="dashed")
    axs[1].set_title(raw"$\chi(t, x)$ (simulação)")

    function ϕ(x)
        if 0 < x < l / 4
            -x^2 / 2
        elseif l / 4 < x < 3l / 4
            x^2 / 2 - x * l / 2 + l^2 / 16
        elseif 3l / 4 < x < l
            -(l - x)^2 / 2
        else
            0.0
        end
    end
    f(x) = ϕ(x - √2 / 2 + l / 2)

    n = 21
    axs[2].plot(x, χ[:, n]; color="black", label="Simulação")
    axs[2].plot(
        x,
        @. (f(x + t[n]) + f(x - t[n])) / 2;
        color="C3",
        linestyle="dashed",
        label="Analítico",
    )
    axs[2].legend(; frameon=true)
    axs[2].set_xlim(0, √2)
    axs[2].set_xlabel("x")
    axs[2].set_title("\$ \\chi($(t[n]), x) \$")

    fig.savefig(plotsdir("toy_model", "quarter-rest.pdf"))
    fig
end

let l = 0.2, V = 0.75, α = 0.0, v₀ = 0.0, model = toy
    x₀ = -(√2 - x_R(α, V; l, v₀) - x_L(α, V; l, v₀)) / 2
    parameters = KinkOscillon(; l, V, α, v₀, x₀, model)
    data, _ = produce_or_load(datadir("toy_model"), parameters) do parameters
        return simulation(parameters; dx=5e-4, sampling=20)
    end
    @unpack x, t, η, H = data
    χ = η - toykink.(t', x)

    χC(t, x) = l^2 / 48

    function χ1R(t, x)
        return (-t + x)^2 +
               (-t + x) * (-√2 + (√7 * l) / 4) +
               (96 - 24 * √14 * l + 25 * l^2) / 192
    end

    function χ2R(t, x)
        return (-5 * (-t + x)^2) / 2 +
               (-t + x) * (5 / √2 - (√7 * l) / 8) +
               (-480 + 24 * √14 * l + 23 * l^2) / 384
    end

    function χ3R(t, x)
        return (-t + x)^2 - (-t + x) * (√2 + (√7 * l) / 4) +
               (96 + 24 * √14 * l + 13 * l^2) / 192
    end

    function χ4R(t, x)
        return (5 * (-t + x)^2) / 2 - (-t + x) * (5 / √2 + (3 * √7 * l) / 8) +
               (160 + 24 * √14 * l + 11 * l^2) / 128
    end

    function χ5R(t, x)
        return -(-t + x)^2 +
               (-t + x) * (√2 + (√7 * l) / 4) +
               (-32 - 8 * √14 * l - 7 * l^2) / 64
    end

    function χ1L(t, x)
        return -1 / 7 * (t + x)^2 +
               ((t + x) * (4 * √2 + √7 * l)) / 28 +
               (-96 - 24 * √14 * l + 7 * l^2) / 1344
    end

    function χ2L(t, x)
        return -1 / 7 * (t + x)^2 +
               ((t + x) * (4 * √2 + √7 * l)) / 28 +
               (-96 - 24 * √14 * l + 7 * l^2) / 1344
    end

    function χ3L(t, x)
        return (t + x)^2 / 14 +
               ((t + x) * (-4 * √2 + √7 * l)) / 56 +
               (32 - 8 * √14 * l + 7 * l^2) / 896
    end

    function χ4L(t, x)
        return (t + x)^2 / 14 +
               ((t + x) * (-4 * √2 + √7 * l)) / 56 +
               (32 - 8 * √14 * l + 7 * l^2) / 896
    end

    function χ5L(t, x)
        return (t + x)^2 / 14 +
               ((t + x) * (-4 * √2 + √7 * l)) / 56 +
               (32 - 8 * √14 * l + 7 * l^2) / 896
    end

    function χsum(t, x)
        xc = √2 / 2 - l / (2γ(V))
        z(x) = xc + l * x
        x0 = 0
        x1 = 2 / (7 * γ(V))
        x2 = 4 / (7 * γ(V))
        x3 = 2 / (3 * γ(V))
        x4 = 6 / (7 * γ(V))
        x5 = 1 / γ(V)

        if z(x5) - t < x ≤ z(x0) + t
            χC(t, x)
        elseif z(x0) + t < x ≤ z(x1) + t
            χ1R(t, x)
        elseif z(x1) + t < x ≤ z(x2) + t
            χ2R(t, x)
        elseif z(x2) + t < x ≤ z(x3) + t
            χ3R(t, x)
        elseif z(x3) + t < x ≤ z(x4) + t
            χ4R(t, x)
        elseif z(x4) + t < x ≤ z(x5) + t
            χ5R(t, x)
        elseif z(x4) - t < x ≤ z(x5) - t
            χ1L(t, x)
        elseif z(x3) - t < x ≤ z(x4) - t
            χ2L(t, x)
        elseif z(x2) - t < x ≤ z(x3) - t
            χ3L(t, x)
        elseif z(x1) - t < x ≤ z(x2) - t
            χ4L(t, x)
        elseif z(x0) - t < x ≤ z(x1) - t
            χ5L(t, x)
        else
            0.0
        end
    end

    fig, axs = plt.subplots(1, 2; figsize=(6.2, 2.8))
    heatmap!(axs[1, 1], x, t, χ; norm=mpl.colors.CenteredNorm(), cmap="RdBu", colorbar=true)
    axs[1].set_xlim(√2 - 5, √2 + 5)
    axs[1].axvline(0; linewidth=0.5, color="black", linestyle="dashed")
    axs[1].axvline(2√2; linewidth=0.5, color="black", linestyle="dashed")
    axs[1].set_title(raw"$\chi(t, x)$ (simulação)")

    n = 21
    axs[2].plot(x, χ[:, n]; color="black", label="Simulação")
    axs[2].plot(x, χsum.(t[n], x); color="C3", linestyle="dashed", label="Analítico")
    axs[2].legend(; frameon=true)
    axs[2].set_xlim(0, √2)
    axs[2].set_xlabel("x")
    axs[2].set_title("\$ \\chi($(t[n]), x) \$")

    fig.savefig(plotsdir("toy_model", "quarter-moving.pdf"))
    fig
end
