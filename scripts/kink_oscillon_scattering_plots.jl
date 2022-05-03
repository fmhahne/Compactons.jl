using Printf
using JLD2, CodecZlib, UnPack
using Compactons

include("../src/plots.jl")
plt.style.use("matplotlibrc")

mkpath("plots/kink_oscillon_scattering")

for l ∈ 0.5:0.5:3.0
    for V ∈ 0.0:1e-1:0.9
        for α ∈ 0.00:0.25:0.75
            basename = @sprintf "kink_oscillon_scattering/l=%.2f,V=%.2f,alpha=%.2f" l V α
            @unpack x, t, η, H, E₁, E₂, E₃ = load("data/$basename.jld2")

            figsize = (3.1, 1.8)
            xlim = (-5, 10)

            # Field
            fig, ax = plt.subplots(figsize=figsize)
            heatmap!(ax, x, t, η)
            ax.set_xlim(xlim...)
            fig.savefig("plots/$basename.pdf")
            fig.clear()

            # Energy density
            fig, ax = plt.subplots(figsize=figsize)
            heatmap!(ax, x, t, H; cmap="magma", norm=mpl.colors.LogNorm(clip=true, vmin=1e-6))
            ax.set_xlim(xlim...)
            fig.savefig("plots/$basename-energy.pdf")
            fig.clear()

            # Perturbation
            χ = η - kink.(t', x)
            fig, ax = plt.subplots(figsize=figsize)
            heatmap!(ax, x, t, χ; norm=mpl.colors.CenteredNorm(), cmap="RdBu")
            ax.set_xlim(xlim...)
            fig.savefig("plots/$basename-perturbation.pdf")
            fig.clear()

            # Energies
            fig, ax = plt.subplots(figsize=figsize)
            ax.plot(t, E₁; label=raw"$E_1$")
            ax.plot(t, E₂; label=raw"$E_2$")
            ax.plot(t, E₃; label=raw"$E_3$")
            ax.legend()
            ax.grid(color="lightgray")
            fig.savefig("plots/$basename-energies.pdf")
            fig.clear()

            plt.close("all")
        end
    end
end
