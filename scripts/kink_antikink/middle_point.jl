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

function plot_middle!(ax, cb, vsave, tsave, η0s)
    middle = hcat(η0s...)
    ax.pcolormesh(vsave, tsave, middle; rasterized=true)

    ax.set_xlabel(raw"$v$")
    ax.set_ylabel(raw"$t$")

    cb.set_array(middle)
    cb.autoscale()

    return nothing
end

dx = 5e-3

fig, axs = plt.subplot_mosaic(
    [["full", "full", "full"], ["A", "B", "C"]]; figsize=(6, 6), layout="constrained"
)
cb = mpl.cm.ScalarMappable()

let vsave = range(; start=0.0, stop=0.5, length=501)
    sampling = 10
    tmax = 50.0
    tsave = 0.0:(sampling * dx):tmax

    η0s = []
    for v in vsave
        η0 = get_sim_middle(; v, dx, tmax, sampling)
        push!(η0s, η0)
    end

    plot_middle!(axs["full"], cb, vsave, tsave, η0s)

    axs["full"].axvline(0.41305; color="red", linewidth=0.5)
    axs["full"].axvline(0.4133; color="red", linewidth=0.5)
    axs["full"].ticklabel_format(; style="plain")
end

let vsave = range(; start=0.41305, stop=0.4133, length=501)
    sampling = 20
    tmax = 100.0
    tsave = 0.0:(sampling * dx):tmax

    η0s = []
    for v in vsave
        η0 = get_sim_middle(; v, dx, tmax, sampling)
        push!(η0s, η0)
    end

    plot_middle!(axs["A"], cb, vsave, tsave, η0s)
    axs["A"].ticklabel_format(; useOffset=false)
    axs["A"].tick_params(; axis="x", labelrotation=90)
    axs["A"].label_outer()

    axs["A"].axvline(0.413245; color="red", linewidth=0.5)
    axs["A"].axvline(0.413250; color="red", linewidth=0.5)
end

let vsave = range(; start=0.413245, stop=0.413250, length=501)
    tmax = 100.0
    sampling = 20
    tsave = 0.0:(sampling * dx):tmax

    η0s = []
    for v in vsave
        η0 = get_sim_middle(; v, dx, tmax, sampling)
        push!(η0s, η0)
    end

    plot_middle!(axs["B"], cb, vsave, tsave, η0s)
    axs["B"].ticklabel_format(; useOffset=false)
    axs["B"].tick_params(; axis="x", labelrotation=90)
    axs["B"].label_outer()

    axs["B"].axvline(0.41324935; color="red", linewidth=0.5)
    axs["B"].axvline(0.41324945; color="red", linewidth=0.5)
end

let vsave = range(; start=0.41324935, stop=0.41324945, length=501)
    tmax = 100.0
    sampling = 20
    tsave = 0.0:(sampling * dx):tmax

    η0s = []
    for v in vsave
        η0 = get_sim_middle(; v, dx, tmax, sampling)
        push!(η0s, η0)
    end

    plot_middle!(axs["C"], cb, vsave, tsave, η0s)
    axs["C"].label_outer()
    axs["C"].ticklabel_format(; useOffset=false)
    axs["C"].tick_params(; axis="x", labelrotation=90)
end

fig.colorbar(cb; ax=[axs["full"], axs["A"], axs["B"], axs["C"]], aspect=40)

Base.mkpath(plotsdir("kink_antikink"))
fig.savefig(plotsdir("kink_antikink", "middle.pdf"))
fig
