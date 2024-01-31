using DrWatson
using Compactons
include(srcdir("plots.jl"))

let
    fig, axs = plt.subplots(
        3,
        2;
        figsize=(6, 6),
        sharex="col",
        sharey="row",
        gridspec_kw=Dict("height_ratios" => [20, 15, 10]),
    )

    for (i, (v, tmax)) in enumerate(zip([0.3, 0.6, 0.9], [20.0, 15.0, 10.0]))
        data, _ = produce_or_load(
            datadir("kink_antikink"), KinkAntikink(; v, tmax), simulation
        )
        @unpack x, t, η, H = data

        heatmap!(axs[i, 1], x, t, η; colorbar=true, norm=mpl.colors.CenteredNorm())
        axs[i, 1].set_xlim(-10, 10)
        axs[i, 1].set_title(raw"$\eta(t,x)\quad v=" * "$v" * raw"$")
        axs[i, 1].label_outer()

        heatmap!(
            axs[i, 2],
            x,
            t,
            H;
            colorbar=true,
            cmap="magma",
            norm=mpl.colors.SymLogNorm(1e-5),
        )
        axs[i, 2].set_xlim(-10, 10)
        axs[i, 2].set_title(raw"$\mathcal{H}(t,x)\quad v=" * "$v" * raw"$")
        axs[i, 2].label_outer()
    end

    Base.mkpath(plotsdir("kink_antikink"))
    fig.savefig(plotsdir("kink_antikink", "examples.pdf"))
    fig
end

let v = 0.3, tmax = 20.0
    fig, axs = plt.subplots(1, 2; figsize=(6, 3), sharex="col", sharey="row")

    data, _ = produce_or_load(datadir("kink_antikink"), KinkAntikink(; v, tmax), simulation)
    @unpack x, t, η, H = data

    heatmap!(axs[1], x, t, η; colorbar=true, norm=mpl.colors.CenteredNorm())
    axs[1].set_xlim(-10, 10)
    axs[1].set_title(raw"$\eta(t,x)$")
    axs[1].label_outer()

    heatmap!(axs[2], x, t, H; colorbar=true, cmap="magma", norm=mpl.colors.SymLogNorm(1e-5))
    axs[2].set_xlim(-10, 10)
    axs[2].set_title(raw"$\mathcal{H}(t,x)$")
    axs[2].label_outer()

    Base.mkpath(plotsdir("kink_antikink"))
    fig.savefig(plotsdir("kink_antikink", "v=$v.pdf"))
    fig
end

let v = 0.6, tmax = 15.0
    fig, axs = plt.subplots(1, 2; figsize=(6, 2.5), sharex="col", sharey="row")

    data, _ = produce_or_load(datadir("kink_antikink"), KinkAntikink(; v, tmax), simulation)
    @unpack x, t, η, H = data

    heatmap!(axs[1], x, t, η; colorbar=true, norm=mpl.colors.CenteredNorm())
    axs[1].set_xlim(-10, 10)
    axs[1].set_title(raw"$\eta(t,x)$")
    axs[1].label_outer()

    heatmap!(axs[2], x, t, H; colorbar=true, cmap="magma", norm=mpl.colors.SymLogNorm(1e-5))
    axs[2].set_xlim(-10, 10)
    axs[2].set_title(raw"$\mathcal{H}(t,x)$")
    axs[2].label_outer()

    Base.mkpath(plotsdir("kink_antikink"))
    fig.savefig(plotsdir("kink_antikink", "v=$v.pdf"))
    fig
end

let v = 0.9, tmax = 10.0
    fig, axs = plt.subplots(1, 2; figsize=(6, 2), sharex="col", sharey="row")

    data, _ = produce_or_load(datadir("kink_antikink"), KinkAntikink(; v, tmax), simulation)
    @unpack x, t, η, H = data

    heatmap!(axs[1], x, t, η; colorbar=true, norm=mpl.colors.CenteredNorm())
    axs[1].set_xlim(-10, 10)
    axs[1].set_title(raw"$\eta(t,x)$")
    axs[1].label_outer()

    heatmap!(axs[2], x, t, H; colorbar=true, cmap="magma", norm=mpl.colors.SymLogNorm(1e-5))
    axs[2].set_xlim(-10, 10)
    axs[2].set_title(raw"$\mathcal{H}(t,x)$")
    axs[2].label_outer()

    Base.mkpath(plotsdir("kink_antikink"))
    fig.savefig(plotsdir("kink_antikink", "v=$v.pdf"))
    fig
end
