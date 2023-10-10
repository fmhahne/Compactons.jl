using DrWatson
using Compactons
include(srcdir("plots.jl"))

let η = -3:1e-3:3, x = -0.5:1e-3:4
    fig, axs = plt.subplots(1, 2; figsize=(6.2, 3))

    for k in 0:0.25:1
        axs[1].plot(η, generalizedmodel(k).V.(η); label="\$k = $k\$")
    end

    axs[1].set_xlabel(raw"$\eta$")
    axs[1].legend(; frameon=true)
    axs[1].set_title(raw"$V_k(\eta)$")
    axs[1].set_xlim(η[begin], η[end])

    for k in 0:0.25:1
        axs[2].plot(x, generalizedkink.(0, x; k=k); label="\$k = $k\$")
    end

    axs[2].set_xlabel(raw"$x$")
    axs[2].set_title(raw"$\eta_K(x)$")
    axs[2].legend(; frameon=true)
    axs[2].set_xlim(x[begin], x[end])

    fig.savefig(plotsdir("generalized_model.pdf"))
    fig
end

let α = 0.25, l = 1.0, V = 0.0, v₀ = 0.0
    fig, ax = plt.subplots()

    for k in 0.98:0.01:1.02
        params = @strdict l V α v₀ k
        data, _ = produce_or_load(datadir("generalized_model"), params) do params
            @unpack l, V, α, v₀, k = params
            model = generalizedmodel(k)

            dx = 5e-4
            x = -3:dx:3

            tsave = 0.0:1e-2:0.75
            tspan = tsave[1], tsave[end]
            xsave = x[begin:10:end]

            η₀ = @. generalizedkink(0, x + x₀(k); k=k) +
                oscillon(α * l, x + l / 2, V; l, v₀)
            ∂ₜη₀ = @. ∂ₜoscillon(α * l, x + l / 2, V; l, v₀)

            η, H = produce_data(model, ∂ₜη₀, η₀, tsave; dx, dt=dx / 10, sampling=10)
            return Dict("x" => xsave, "t" => tsave, "η" => η, "H" => H)
        end

        @unpack x, t, η, H = data
        χ = η - generalizedkink.(t', x .+ x₀(k); k=k)

        ax.plot(x, χ[:, end] - χ[end:-1:1, end]; label="\$k=$k\$")
    end

    ax.set_xlim(-1.5, 1.5)
    ax.set_title(raw"$\chi(0.75, x) - \chi(0.75, -x)$")
    ax.set_xlabel(raw"$x$")
    ax.legend(; frameon=true)

    fig.savefig(plotsdir("assymmetry.pdf"))
    fig
end
