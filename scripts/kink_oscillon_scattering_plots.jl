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
            @unpack t, x, η, H = load("data/" * basename * ".jld2")

            figsize = (3.1, 2.6)
            xlim = (-t[end], t[end]) ./ 2

            # Field
            fig, ax = plt.subplots(figsize=figsize)
            heatmap!(ax, x, t, η)
            ax.set_xlim(xlim...)
            fig.savefig("plots/" * basename * ".pdf")
            fig.clear()

            # Energy density
            fig, ax = plt.subplots(figsize=figsize)
            heatmap!(ax, x, t, H; cmap="magma", norm=mpl.colors.LogNorm(clip=true, vmin=1e-6))
            ax.set_xlim(xlim...)
            fig.savefig("plots/" * basename * "-energy.pdf")
            fig.clear()

            # Perturbation
            χ = η - kink.(t', x .+ π / γ(V), V)
            fig, ax = plt.subplots(figsize=figsize)
            heatmap!(ax, x, t, χ; norm=mpl.colors.CenteredNorm(), cmap="RdBu")
            ax.set_xlim(xlim...)
            fig.savefig("plots/" * basename * "-perturbation.pdf")
            fig.clear()

            plt.close("all")
        end
    end
end
