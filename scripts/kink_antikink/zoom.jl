using Compactons
using DrWatson
include(srcdir("plots.jl"))

let v = 0.75
    fig, axs = plt.subplots(1, 2; figsize=(6, 2.78))

    params = KinkAntikink(; v)
    data, _ = produce_or_load(datadir("kink_antikink"), params, simulation)
    @unpack x, t, η, H = data

    x1 = -1.5
    x2 = 1.5
    t1 = 2.5
    t2 = 5.5

    heatmap!(axs[1], x, t, H; colorbar=true, cmap="magma", norm=mpl.colors.SymLogNorm(1e-5))
    axs[1].set_xlim(-5, 5)
    axs[1].add_patch(
        mpl.patches.Rectangle(
            (x1, t1), x2 - x1, t2 - t1; facecolor="none", edgecolor="lime", linewidth=0.5
        ),
    )
    axs[1].set_title(raw"$\mathcal{H}(t, x)$")

    nx1 = findfirst(x .≥ x1)
    nx2 = findfirst(x .≥ x2)
    nt1 = findfirst(t .≥ t1)
    nt2 = findfirst(t .≥ t2)

    vlim = maximum(abs.(η[nx1:nx2, nt1:nt2] .+ 2.0))
    heatmap!(
        axs[2],
        x[nx1:nx2],
        t[nt1:nt2],
        η[nx1:nx2, nt1:nt2] .+ 2.0;
        colorbar=true,
        cmap="RdBu",
        norm=mpl.colors.SymLogNorm(1e-3; vmin=-vlim, vmax=vlim),
    )
    axs[2].set_title(raw"$\eta(t, x) + 2$")

    Base.mkpath(plotsdir("kink_antikink"))
    fig.savefig(plotsdir("kink_antikink", "zoom.pdf"))
    fig
end
