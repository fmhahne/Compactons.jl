using DrWatson
using Compactons
include(srcdir("plots.jl"))

vsave = 0.0:0.001:0.5
sampling = 10
dx = 5e-3
saveat = dx * sampling
tmax = 50.0
tsave = 0.0:saveat:tmax

let η0s = []
    for v in vsave
        params = KinkAntikinkMiddle(; v, dx, tmax)
        data, _ = produce_or_load(datadir("kink_antikink", "middle"), params, simulation)
        @unpack η0 = data
        push!(η0s, η0)
    end

    middle = hcat(η0s...)

    fig, ax = plt.subplots()
    img = ax.pcolormesh(vsave, tsave, middle; rasterized=true)
    ax.set_xlabel(raw"$v$")
    ax.set_ylabel(raw"$t$")
    fig.colorbar(img)

    Base.mkpath(plotsdir("kink_antikink"))
    fig.savefig(plotsdir("kink_antikink", "middle.pdf"))
    fig
end

let η0s = []
    middle = zero(tsave * vsave')
    for (i, v) in collect(enumerate(vsave))
        params = CCKinkAntikink(; v, tmax, saveat)
        data, _ = produce_or_load(
            datadir("cc_kink_antikink"), params, collective_coordinates
        )
        @unpack a, c = data
        push!(η0s, ηKAK.(0.0, a, c))
    end

    middle = hcat(η0s...)

    fig, ax = plt.subplots()
    img = ax.pcolormesh(vsave, tsave, middle; rasterized=true)
    ax.set_xlabel(raw"$v$")
    ax.set_ylabel(raw"$t$")
    fig.colorbar(img)

    Base.mkpath(plotsdir("kink_antikink"))
    fig.savefig(plotsdir("kink_antikink", "middle_point_cc.pdf"))
    fig
end
