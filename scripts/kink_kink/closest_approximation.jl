using DrWatson
using Compactons
include(srcdir("plots.jl"))

tmax = 20.0
vsave0 = 0.0:0.002:0.999
vsave1 = 0.0:0.002:0.70
vsave2 = 0.0:0.002:0.87
xRmins = Float64[]
xRmins_non_rel = Float64[]
xRmins_non_rel_mode = Float64[]
xRmins_rel = Float64[]
xRmins_rel_mode = Float64[]

for v in vsave0
    params = KinkKinkBorder(; v, xmax=tmax, tmax)
    data, _ = produce_or_load(datadir("kink_kink_border"), params, simulation)
    @unpack x, xR_idxs = data
    xRs = Float64[]
    for i in xR_idxs.saveval
        append!(xRs, x[i])
    end
    append!(xRmins, minimum(xRs))
end

for v in vsave1
    params = CCKinkKinkNonRel(; v, tmax)
    data, _ = produce_or_load(
        datadir("cc_kink_kink_non_rel"), params, collective_coordinates
    )
    @unpack a = data
    append!(xRmins_non_rel, π / 2 + minimum(a))

    params = CCKinkKinkNonRelMode(; v, tmax)
    data, _ = produce_or_load(
        datadir("cc_kink_kink_non_rel_mode"), params, collective_coordinates
    )
    @unpack a = data
    append!(xRmins_non_rel_mode, π / 2 + minimum(a))
end

for v in vsave2
    params = CCKinkKinkRel(; v, tmax)
    data, _ = produce_or_load(datadir("cc_kink_kink_rel"), params, collective_coordinates)
    @unpack a, b = data
    append!(xRmins_rel, minimum(@. π / (2 * b) + a))

    params = CCKinkKinkRelMode(; v, tmax)
    data, _ = produce_or_load(
        datadir("cc_kink_kink_rel_mode"), params, collective_coordinates
    )
    @unpack a, b = data
    append!(xRmins_rel_mode, minimum(@. π / (2 * b) + a))
end

fig, ax = plt.subplots(; figsize=(4.5, 3.0))
ax.plot(vsave0, xRmins; label="Simulation")
ax.plot(vsave1, xRmins_non_rel; label="Non-relativistic")
ax.plot(vsave1, xRmins_non_rel_mode; label="Non-relativistic with internal mode")
ax.plot(vsave2, xRmins_rel; label="Relativistic")
ax.plot(vsave2, xRmins_rel_mode; label="Relativistic with internal mode")
ax.set_xlabel(raw"$v$")
ax.set_ylabel(raw"$x_{R,\mathrm{min}}$")
ax.legend(; frameon=true)
fig.savefig(plotsdir("kink_kink", "closest_approximation.pdf"))
fig
