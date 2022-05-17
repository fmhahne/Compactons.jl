using DrWatson
using Compactons
include(srcdir("plots.jl"))

let η = -3:1e-3:3, x = -0.5:1e-3:4
    fig, axs = plt.subplots(1, 2; figsize=(6.2, 3))

    for k ∈ 0:0.25:1
        axs[1].plot(η, generalizedmodel(k).V.(η); label="\$k = $k\$")
    end

    axs[1].set_xlabel(raw"$\eta$")
    axs[1].legend(frameon=true)
    axs[1].set_title(raw"$V_k(\eta)$")
    axs[1].set_xlim(η[begin], η[end])

    for k ∈ 0:0.25:1
        axs[2].plot(x, generalizedkink.(0, x; k=k); label="\$k = $k\$")
    end

    axs[2].set_xlabel(raw"$x$")
    axs[2].set_title(raw"$\eta_K(x)$")
    axs[2].legend(frameon=true)
    axs[2].set_xlim(x[begin], x[end])

    fig.savefig(plotsdir("generalized_model.pdf"))
    fig
end

let
    dx = 5e-4
    x = -3:dx:3
    N = length(x)

    tsave = 0.0:1e-2:0.75
    tspan = tsave[1], tsave[end]

    xsave = x[1:10:N]

    l = 1
    Δ = 0
    V = 0

    fig, ax = plt.subplots()
    for k ∈ 0.98:0.01:1.02
        model = generalizedmodel(k)
        η₀ = generalizedkink.(0, x .+ x₀(k); k=k) + oscillon.(0.25 * l, x .+ 0.5 * l .- Δ, V; l=l)
        ∂ₜη₀ = ∂ₜoscillon.(0.25 * l, x .+ 0.5 * l .- Δ, V; l=l)

        η, H = producedata(model, ∂ₜη₀, η₀, tsave; dx, dt=dx / 10, sampling=10)
        χ = η - generalizedkink.(tsave', xsave .+ x₀(k); k=k)
        ax.plot(xsave, χ[:, end] - χ[end:-1:1, end]; label="\$k=$k\$")
    end
    ax.set_xlim(-1.5, 1.5)
    ax.set_title(raw"$\chi(0.75, x) - \chi(0.75, -x)$")
    ax.set_xlabel(raw"$x$")
    ax.legend(frameon=true)

    fig.savefig(plotsdir("assymmetry.pdf"))
    fig
end
