using Compactons
using DrWatson
include(srcdir("plots.jl"))

a = (-π / 2):0.01:(π / 2)
c = -1:0.01:1

U = @. potential(
    a',
    c;
    model=quadratic,
    η=ηKK_non_rel_mode,
    xspan=xspanKK_non_rel,
    quadgk_kwargs=(; atol=1e-6, order=20),
)

fig, axs = plt.subplots(1, 2; figsize=(6, 2.78))

img = axs[1].contourf(a, c, U)
plt.colorbar(img)
axs[1].set_xlabel(raw"$a$")
axs[1].set_ylabel(raw"$c$")
axs[1].set_xlim(a[begin], a[end])
axs[1].set_ylim(c[begin], c[end])

for (v, color, linestyle) in zip([0.557, 0.556], ["red", "lime"], ["dashed", "solid"])
    params = CCKinkKinkNonRelMode(; v, tmax=20.0)
    data, _ = produce_or_load(
        datadir("cc_kink_kink_non_rel_mode"), params, collective_coordinates
    )
    axs[1].plot(data["a"], data["c"]; color, linestyle, label="\$v=$v\$")
end
axs[1].legend(; frameon=true)

for i in [101, 96, 91]
    axs[2].plot(a, U[i, :]; label="\$U(a, $(c[i]))\$")
end
axs[2].legend()
axs[2].set_xlabel(raw"$a$")

mkpath(plotsdir("kink_kink"))
fig.savefig(plotsdir("kink_kink", "non_rel_mode_potential.pdf"))
fig
