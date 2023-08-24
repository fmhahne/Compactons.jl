using DrWatson
using FFTW
using Compactons
include(srcdir("plots.jl"))

function C(t, k)
    if k^2 ≥ 1
        cos(t * √(k^2 - 1))
    else
        cosh(t * √(1 - k^2))
    end
end

function S(t, k)
    if k^2 ≥ 1
        sin(t * √(k^2 - 1)) / √(k^2 - 1)
    else
        sinh(t * √(1 - k^2)) / √(1 - k^2)
    end
end

"""
    χ(t, k, f̃ₖ, g̃ₖ)

Perturbation calculated from the Fourier transforms `f̃ₖ` and `g̃ₖ` of the initial conditions,
using frequencies `k`.
"""
χ(t, k, f̃ₖ, g̃ₖ) = ifft(f̃ₖ .* C.(t, k) .+ g̃ₖ .* S.(t, k))

let l = 1.0, α = 0.75, v₀ = 0.0, ns = [31, 61, 91]
    fig, axs = plt.subplots(2, length(ns); figsize=(6.2, 6.2 / 1.62), tight_layout=false, sharex="col", sharey="row")

    for (i, V) in enumerate([0, 0.75])
        x₀ = -(π - x_R(α, V; l, v₀) - x_L(α, V; l, v₀)) / 2
        data, _ = produce_or_load(datadir("kink_oscillon_superposition"), KinkOscillon(; l, V, α, v₀, x₀), simulation)
        @unpack t, x, η, H = data

        k = fftfreq(length(x)) * length(x) / π
        x₀ = -(π - x_R(α, V; l, v₀) - x_L(α, V; l, v₀)) / 2
        f̃ₖ = fft(@. oscillon.(l * α * γ(V), x + x₀, V; l, v₀))
        g̃ₖ = fft(@. ∂ₜoscillon(l * α * γ(V), x + x₀, V; l, v₀))

        for (n, ax) in zip(ns, axs[i, 1:3])
            ax.set_title("\$V = $V, \\, t = $(t[n])\$")
            ax.plot(x, η[:, n] - kink.(x); label="Simulação", color="black")
            ax.plot(x, real.(χ(t[n], k, f̃ₖ, g̃ₖ)); label="Semi-analítico", color="C3", linestyle="dashed")
            ax.set_xlim(0, float(π))
            ax.set_xlabel(raw"$x$")
            ax.label_outer()
        end
    end

    fig.tight_layout()

    handles, labels = axs[1, 1].get_legend_handles_labels()
    fig.subplots_adjust(bottom=0.18)
    fig.legend(handles, labels; loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0))

    fig.savefig(plotsdir("kink_oscillon_superposition", "comparison.pdf"))
    fig
end
