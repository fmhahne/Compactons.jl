using DrWatson
using Compactons
include(srcdir("plots.jl"))

vsave = 0.0:0.001:0.499
sampling = 10
dx = 5e-3
saveat = dx * sampling
tmax = 50.0
tsave = 0.0:saveat:tmax

η0s = []
for v in vsave
    data, _ = produce_or_load(
        datadir("kink_antikink"), KinkAntikink(; v, dx, tmax, sampling), simulation
    )
    @unpack η = data
    push!(η0s, η[end ÷ 2 + 1, :])
end

middle = hcat(η0s...)

fig, ax = plt.subplots()
img = ax.pcolormesh(vsave, tsave, middle; rasterized=true)
ax.set_xlabel(raw"$v$")
ax.set_ylabel(raw"$t$")
fig.colorbar(img)

Base.mkpath(plotsdir("kink_antikink"))
fig.savefig(plotsdir("kink_antikink", "middle_point.pdf"))
fig
