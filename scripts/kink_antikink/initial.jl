using DrWatson
using Compactons
include(srcdir("plots.jl"))
mkpath(plotsdir("kink_antikink"))

let v = 0.5
    fig, ax = plt.subplots()
    x = range(-π, π, 1000)

    η₀ = kink.(0.0, -abs.(x) .+ π / γ(v), v)
    ∂ₜη₀ = ∂ₜkink.(0, -abs.(x) .+ π / γ(v), v)

    ax.plot(x, η₀; color="black", label=raw"$\eta(0, x)$")
    ax.plot(x, ∂ₜη₀; color="gray", ls="dashed", label=raw"$\partial_t \eta(0, x)$")
    ax.legend()
    fig.savefig(plotsdir("kink_antikink", "initial.pdf"))
    fig
end
