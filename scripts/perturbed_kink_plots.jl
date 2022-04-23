using Printf
using JLD2, CodecZlib, UnPack

include("../src/plots.jl")
plt.style.use("matplotlibrc")

mkpath("plots/perturbed_kink")

for ϵ ∈ -0.50:0.05:0.50
    basename = @sprintf "perturbed_kink/eps=%.2f" ϵ
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

    plt.close("all")
end
