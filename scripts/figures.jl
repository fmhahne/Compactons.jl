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

    for ((α, v₀), ax) in zip(Iterators.product(αs, v₀s), axs)
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

let l = 1.0, v₀ = 0.6, Vs = [0.0, 0.4]
    fig, axs = plt.subplots(1, 2; figsize=(6.2, 3.1), tight_layout=false, sharey=true)
    norm = mpl.colors.CenteredNorm()
    cmap = "RdBu"
    cb = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

    Δ = v₀ * l / 2
    for (ax, V) in zip(axs, Vs)
        x = (-V*l):1e-3:(l+Δ+V*l)
        t = -l:1e-3:l
        φ = oscillon.(t', x, V; l, v₀)

        heatmap!(ax, x, t, φ; cmap=cmap, norm=norm)
        ax.set_title("\$ V = $V \$")
        ax.label_outer()

        cb.set_array(φ)
        cb.autoscale()
    end

    shared_colorbar!(fig, cb)
    fig.savefig(plotsdir("oscillon.pdf"))
    fig
end

let x = -0.1:1e-3:(π+0.1)
    fig, ax = plt.subplots()
    ax.plot(x, kink.(x); label=raw"$\eta_K(x)$", color="black")
    ax.plot(x, @. 2 - kink(x); label=raw"$\eta_\bar{K}(x)$", color="gray", linestyle="dashed")
    ax.legend()
    ax.set_xlim(x[begin], x[end])
    fig.savefig(plotsdir("kink.pdf"))
    fig
end
