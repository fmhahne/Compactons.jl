using Printf
using JLD2, CodecZlib, UnPack

include("../src/plots.jl")
plt.style.use("matplotlibrc")

mkpath("plots/kink_antikink_scattering")

for V ∈ 0.00:1e-2:0.99
    basename = @sprintf "kink_antikink_scattering/V=%.2f" V
    @unpack t, x, η, H = load("data/$basename.jld2")

    figsize = V < 0.5 ? (3.1, 2.6) : (3.1, 1.55)
    xlim = -10, 10

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

    plt.close("all")
end
