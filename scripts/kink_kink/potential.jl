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

fig, ax = plt.subplots()
img = ax.contourf(a, c, U)
ax.set_xlabel(raw"$a$")
ax.set_ylabel(raw"$c$")
plt.colorbar(img)
fig.savefig(plotsdir("kink_kink", "non_rel_mode_potential.pdf"))
fig
