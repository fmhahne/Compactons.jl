using DrWatson
using Compactons
include(srcdir("plots.jl"))

tmax = 20.0
vsave = 0:1e-4:0.01
xRmins = Float64[]
xRmins_non_rel = Float64[]
xRmins_non_rel_mode = Float64[]
xRmins_rel = Float64[]
xRmins_rel_mode = Float64[]

for v in vsave
    params = KinkKinkBorder(; v, xmax=5.0, tmax)
    data, _ = produce_or_load(datadir("kink_kink_border"), params, simulation)
    @unpack x, xR_idxs = data

    xRs = Float64[]
    for i in xR_idxs.saveval
        append!(xRs, x[i])
    end
    append!(xRmins, minimum(xRs))

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

fig, ax = plt.subplots()
ax.scatter(vsave, xRmins; color="black", marker="o", s=3, label="Simulation")
ax.scatter(vsave, xRmins_non_rel; marker="x", s=3, label="Non-relativistic")
ax.scatter(
    vsave, xRmins_non_rel_mode; marker="+", s=3, label="Non-relativistic with internal mode"
)
ax.scatter(vsave, xRmins_rel; marker="D", s=3, label="Relativistic")
ax.scatter(vsave, xRmins_rel_mode; marker="^", s=3, label="Relativistic with internal mode")
ax.set_xlabel(raw"$v$")
ax.set_ylabel(raw"$x_{R,\mathrm{min}}$")
ax.legend(; frameon=true)
fig.savefig(plotsdir("kink_kink", "closest_approximation.pdf"))
fig
