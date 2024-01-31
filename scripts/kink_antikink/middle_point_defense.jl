using DrWatson
using Compactons
include(srcdir("plots.jl"))

function get_sim_middle(; v, dx, tmax, sampling)
    params = KinkAntikinkMiddle(; v, dx, tmax, sampling)
    data, _ = produce_or_load(
        datadir("kink_antikink", "middle"),
        params,
        simulation;
        filename=params -> savename(params; sigdigits=12),
    )
    return data["η0"]
end

function plot_middle!(fig, ax, vsave, tsave, η0s)
    middle = hcat(η0s...)
    im = ax.pcolormesh(vsave, tsave, middle; rasterized=true)
    fig.colorbar(im; ax)

    ax.set_xlabel(raw"$v$")
    ax.set_ylabel(raw"$t$")
    return nothing
end

dx = 5e-3

let vsave = range(; start=0.0, stop=0.5, length=501)
    fig, ax = plt.subplots()

    sampling = 10
    tmax = 50.0
    tsave = 0.0:(sampling * dx):tmax

    η0s = []
    for v in vsave
        η0 = get_sim_middle(; v, dx, tmax, sampling)
        push!(η0s, η0)
    end

    plot_middle!(fig, ax, vsave, tsave, η0s)

    ax.axvline(0.41305; color="red", linewidth=0.5)
    ax.axvline(0.4133; color="red", linewidth=0.5)

    fig.savefig(plotsdir("kink_antikink", "middle0.pdf"))
    fig
end

let vsave = range(; start=0.41305, stop=0.4133, length=501)
    fig, ax = plt.subplots()

    sampling = 20
    tmax = 100.0
    tsave = 0.0:(sampling * dx):tmax

    η0s = []
    for v in vsave
        η0 = get_sim_middle(; v, dx, tmax, sampling)
        push!(η0s, η0)
    end

    plot_middle!(fig, ax, vsave, tsave, η0s)

    ax.axvline(0.413245; color="red", linewidth=0.5)
    ax.axvline(0.413250; color="red", linewidth=0.5)

    fig.savefig(plotsdir("kink_antikink", "middle1.pdf"))
    fig
end

let vsave = range(; start=0.413245, stop=0.413250, length=501)
    fig, ax = plt.subplots()

    tmax = 100.0
    sampling = 20
    tsave = 0.0:(sampling * dx):tmax

    η0s = []
    for v in vsave
        η0 = get_sim_middle(; v, dx, tmax, sampling)
        push!(η0s, η0)
    end

    plot_middle!(fig, ax, vsave, tsave, η0s)

    ax.axvline(0.41324935; color="red", linewidth=0.5)
    ax.axvline(0.41324945; color="red", linewidth=0.5)

    ax.ticklabel_format(; useOffset=false)
    ax.tick_params(; axis="x", labelrotation=90)

    fig.savefig(plotsdir("kink_antikink", "middle2.pdf"))
    fig
end

let vsave = range(; start=0.41324935, stop=0.41324945, length=501)
    fig, ax = plt.subplots()

    tmax = 100.0
    sampling = 20
    tsave = 0.0:(sampling * dx):tmax

    η0s = []
    for v in vsave
        η0 = get_sim_middle(; v, dx, tmax, sampling)
        push!(η0s, η0)
    end

    plot_middle!(fig, ax, vsave, tsave, η0s)

    ax.ticklabel_format(; useOffset=false)
    ax.tick_params(; axis="x", labelrotation=90)

    fig.savefig(plotsdir("kink_antikink", "middle3.pdf"))
    fig
end
