using DrWatson
using Compactons
include(srcdir("plots.jl"))

vsave = 0.0:0.001:0.4
sampling = 10
dx = 5e-3
saveat = dx * sampling
tmax = 50.0
tsave = 0.0:saveat:tmax

let Es = []
    for v in 0:1e-3:0.4
        data, _ = produce_or_load(
            datadir("kink_antikink"), KinkAntikink(; v, dx, tmax, sampling), simulation
        )
        @unpack H, x = data
        i = findfirst(x .> π)
        Δx = (x[end] - x[begin]) / (length(x) - 1)
        E = vec(2 * Δx * sum(H[i:end, :]; dims=1)) / γ(v) / π
        push!(Es, E)
    end

    Es = hcat(Es...)

    fig, ax = plt.subplots()
    img = ax.pcolormesh(collect(0:1e-3:0.4), tsave, Es; rasterized=true, cmap="magma")
    ax.set_xlabel(raw"$v$")
    ax.set_ylabel(raw"$t$")
    fig.colorbar(img)

    Base.mkpath(plotsdir("kink_antikink"))
    fig.savefig(plotsdir("kink_antikink", "radiation_energy.pdf"))
    fig
end
