function moduli_space_hamiltonian(p, q, (metric, potential), t)
    return dot(p, inv(metric(q)), p) / 2 + potential(q)
end

function moduli_space_metric(q; model, η, xspan, quadgk_kwargs=(;))
    ∂η_∂q(x) = FiniteDiff.finite_difference_gradient(q -> η(x, q), q)
    xi, xf = xspan(q)
    g, _ = quadgk(x -> ∂η_∂q(x) * ∂η_∂q(x)', xi, xf; quadgk_kwargs...)
    return g
end

function moduli_space_potential(q; model, η, xspan, quadgk_kwargs=(;))
    ∂η_∂x(x) = FiniteDiff.finite_difference_derivative(x -> η(x, q), x)
    xi, xf = xspan(q)
    U, _ = quadgk(x -> abs2(∂η_∂x(x)) + model.V(η(x, q)), xi, xf; quadgk_kwargs...)
    return U
end
