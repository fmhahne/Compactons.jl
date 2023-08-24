using DrWatson
using DifferentialEquations
using Compactons
include(srcdir("plots.jl"))

Vsave = 0.0:0.001:0.499
dx = 5e-3
sampling = 10
tmax = 50.0
tsave = 0.0:(dx * sampling):tmax

η₀s = []
for V in Vsave
    data, _ = produce_or_load(
        datadir("kink_antikink"), KinkAntikink(; V, dx, tmax, sampling), simulation
    )
    @unpack η = data
    push!(η₀s, η[end ÷ 2 + 1, :])
end

data = reduce(hcat, η₀s)

fig, ax = plt.subplots()
img = ax.pcolormesh(Vsave, tsave, data; rasterized=true)
ax.set_xlabel(raw"$V$")
ax.set_ylabel(raw"$t$")
fig.colorbar(img)
fig.savefig(plotsdir("kink_antikink", "middle-point.pdf"))
fig
