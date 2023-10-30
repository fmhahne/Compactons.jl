using Compactons
using DrWatson
include(srcdir("plots.jl"))

η0s = []
vsave = range(; start=0.0, stop=0.5, length=501)
saveat = 0.05
tmax = 50.0
tsave = 0.0:saveat:tmax

for (i, v) in collect(enumerate(vsave))
    params = CCKinkAntikink(; v, tmax, saveat)
    data, _ = produce_or_load(datadir("cc_kink_antikink"), params, collective_coordinates)
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
