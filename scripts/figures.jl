using DrWatson
using Compactons
include(srcdir("plots.jl"))

let η = -3:1e-3:3
    fig, ax = plt.subplots()

    ax.plot(η, quadratic.V.(η); label=raw"$V(\eta)$", color="black")
    ax.plot(η, quadratic.V′.(η); label=raw"$V'(\eta)$", color="gray", linestyle="dashed")

    ax.set_xlabel(raw"$\eta$")
    ax.legend()

    fig.savefig(plotsdir("quadratic-potential.pdf"))
    fig
end

let η = -3:1e-3:3
    fig, ax = plt.subplots()

    ax.plot(η, toy.V.(η); label=raw"$V(\eta)$", color="black")
    ax.plot(η, toy.V′.(η); label=raw"$V'(\eta)$", color="gray", linestyle="dashed")

    ax.set_xlabel(raw"$\eta$")
    ax.legend()

    fig.savefig(plotsdir("toy-potential.pdf"))
    fig
end

let l = 2, V = 0.6, v₀s = [0, 0.5], αs = [0, 0.25]
    x = -5:1e-3:10
    fig, axs = plt.subplots(2, 2; figsize=(6.2, 6.2 * 2 / (1 + √5)), sharex=true, sharey=true, tight_layout=true)

    for ((α, v₀), ax) ∈ zip(Iterators.product(αs, v₀s), axs)
        η₀ = kink.(0.0, x) + oscillon.(l * α * γ(V), x .+ x_R(α, V; l, v₀), V; l, v₀)
        ∂ₜη₀ = ∂ₜkink.(0.0, x) + ∂ₜoscillon.(l * α * γ(V), x .+ x_R(α, V; l, v₀), V; l, v₀)

        ax.set_title("\$v_0 = $(v₀)\$, \$\\alpha = $α\$")
        ax.plot(x, η₀; color="black", label=raw"$\eta(0, x)$")
        ax.plot(x, ∂ₜη₀; color="gray", ls="dashed", label=raw"$\partial_t \eta(0, x)$")
        ax.set_xlabel(raw"$x$")
        ax.label_outer()
        ax.set_xlim(-3, 4)
        ax.legend()
    end

    fig.savefig(plotsdir("kink_oscillon_scattering", "initial.pdf"))
    fig
end
