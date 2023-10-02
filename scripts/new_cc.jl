using Compactons
using DifferentialEquations
using DrWatson
include(srcdir("plots.jl"))

dt = 1e-2
tmax = 50.0
ts = 0:dt:tmax
vs = 0.001:0.001:0.5

middle = zeros(Float64, (length(ts), length(vs)))
for (i, v) in enumerate(vs)
    params = KinkAntikinkNonRelModuliSpace(; v, dt, tmax)
    data, _ = produce_or_load(
        datadir("moduli_space", "kink_antikink_non_relativistic"), params, moduli_space
    )
    @unpack a = data

    middle[:, i] = @.(kink(a - π / 2) - kink(-a + π / 2))
end

fig, ax = plt.subplots()
c = ax.pcolormesh(vs, ts, middle; rasterized=true)
fig.colorbar(c; ax=ax)
ax.set_xlabel(raw"$v$")
ax.set_ylabel(raw"$t$")

Base.mkpath(plotsdir("kink_antikink"))
fig.savefig(plotsdir("kink_antikink", "nonrel_middle.pdf"))
fig
