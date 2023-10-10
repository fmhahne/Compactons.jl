using DrWatson
using Compactons
include(srcdir("plots.jl"))

tmax = 20.0
vsave = 0:1e-4:0.01
xRmins = Float64[]
xRmins_a = Float64[]
xRmins_ab = Float64[]

for v in vsave
    params = KinkKinkBorder(; v, xmax=5.0, tmax)
    data, _ = produce_or_load(datadir("kink_kink_border"), params, simulation)
    @unpack x, xR_idxs = data

    xRs = Float64[]
    for i in xR_idxs.saveval
        append!(xRs, x[i])
    end
    append!(xRmins, minimum(xRs))

    cca_params = KKa(; v, tmax)
    data, _ = produce_or_load(datadir("KKa"), cca_params, collectivecoordinates)
    @unpack solution = data
    append!(xRmins_a, π / 2 + minimum(solution[2, :]))

    ccab_params = KKab(; v, tmax)
    data, _ = produce_or_load(datadir("KKab"), ccab_params, collectivecoordinates)
    @unpack solution = data
    append!(xRmins_ab, minimum(@. π / 2 / solution[4, :] + solution[3, :]))
end

fig, ax = plt.subplots()
ax.scatter(vsave, xRmins; color="black", marker="o", s=3, label="Simulation")
ax.scatter(vsave, xRmins_a; marker="+", s=3, label="Non-relativistic")
ax.scatter(vsave, xRmins_ab; marker="x", s=3, label="Relativistic")
ax.set_xlabel(raw"$v$")
ax.set_ylabel(raw"$x_{R,\mathrm{min}}$")
ax.legend()
fig.savefig(plotsdir("kink_kink", "closest_approximation.pdf"))
fig
