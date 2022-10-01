using DrWatson
using DifferentialEquations
using Roots
include(srcdir("plots.jl"))

function vortexeq(u, p, r)
    f′, α′, f, α = u
    λ, e, n = p

    f″ = -f′ / r + (n - α)^2 * f / r^2 - λ * f
    α″ = α′ / r - 2 * e^2 * (n - α) * f^2

    return [f″, α″, f′, α′]
end

function bc!(residual, u, p, r)
    _, _, n = p

    zeros = find_zeros(r -> u(r)[1], r[begin], r[end])

    if length(zeros) == 0
        R = r[end]
    else
        R = zeros[begin]
    end

    residual[1] = u(0.0)[3]
    residual[2] = u(0.0)[2]
    residual[3] = u(R)[3] - 1.0
    residual[4] = u(R)[4] - 1.0

    return nothing
end

λ = 1.0
e = 0.1
n = 1.0
p = λ, e, n

r_max = 5.0
ϵ = 1e-3

f1 = 1.167
a2 = 0.5041

f3 = -f1 * (2a2 + λ) / 8
f5 = f1 * (12 * a2^2 + 4 * e^2 * f1^2 + 4 * a2 * λ + λ^2) / 192

a4 = -e^2 * f1^2 / 4
a6 = e^2.0f1^2 * (6a2 + λ) / 48

fϵ = f1 * ϵ + f3 * ϵ^3 + f5 * ϵ^5
f′ϵ = f1 + 3 * f3 * ϵ^2 + 5 * f5 * ϵ^4

αϵ = a2 * ϵ^2 + a4 * ϵ^4 + a6 * ϵ^6
α′ϵ = 2 * a2 * ϵ + 4 * a4 * ϵ^3 + 6 * a6 * ϵ^5

u₀ = [f′ϵ, α′ϵ, fϵ, αϵ]

prob = BVProblem(vortexeq, bc!, u₀, (ϵ, r_max), p)
sol = solve(prob)

R = find_zeros(r -> sol(r)[1], ϵ, r_max)[begin]
r = 0:1e-2:r_max
f = reduce(hcat, sol.(r))[3, :]
α = reduce(hcat, sol.(r))[4, :]

fig, ax = plt.subplots()

ax.plot(r, f; label=raw"$f\,(r)$")
ax.plot(r, α; label=raw"$\alpha(r)$")
ax.set_xlabel(raw"$r$")
ax.set_ylim(-0.1, 1.1)
ax.set_xlim(0, R)
ax.legend()

fig.savefig(plotsdir("vortex.pdf"))
fig
