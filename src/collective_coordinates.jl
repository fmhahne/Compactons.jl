function cceq!(q̈, q̇, q, (model, η, xspan, quadgk_kwargs), t)
    ∂η_∂q(x) = FiniteDiff.finite_difference_gradient(q -> η(x, q), q)
    ∂η_∂x(x) = FiniteDiff.finite_difference_derivative(x -> η(x, q), x)
    ∂²η_∂q∂q(x) = FiniteDiff.finite_difference_hessian(q -> η(x, q), q)
    ∂²η_∂q∂x(x) = FiniteDiff.finite_difference_derivative(∂η_∂q, x)

    xi, xf = xspan(q)
    g, _ = quadgk(x -> ∂η_∂q(x) * ∂η_∂q(x)', xi, xf; quadgk_kwargs...)
    f, _ = quadgk(
        x ->
            ∂²η_∂q∂x(x) * ∂η_∂x(x) +
            (model.V′(η(x, q)) + dot(q̇', ∂²η_∂q∂q(x), q̇)) * ∂η_∂q(x),
        xi,
        xf;
        quadgk_kwargs...,
    )
    q̈ .= -g \ f
    return nothing
end

function collective_coordinates(params)
    return error("Collective coordinates $(typeof(params)) not implemented")
end

function metric(q...; η, xspan, quadgk_kwargs)
    ∂η_∂q(x) = FiniteDiff.finite_difference_gradient(q -> η(x, q), q)

    xi, xf = xspan(q)
    g, _ = quadgk(x -> ∂η_∂q(x) * ∂η_∂q(x)', xi, xf; quadgk_kwargs...)
    return g
end

function potential(q...; model, η, xspan, quadgk_kwargs)
    ∂η_∂x(x) = FiniteDiff.finite_difference_derivative(x -> η(x, q), x)

    xi, xf = xspan(q)
    U, _ = quadgk(
        x -> 0.5 * dot(∂η_∂x(x), ∂η_∂x(x)) + model.V(η(x, q)), xi, xf; quadgk_kwargs...
    )
    return U
end
