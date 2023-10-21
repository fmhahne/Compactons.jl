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

# Kink-kink scattering

# Position mode

@with_kw struct KKa{T<:Real}
    v::T
    tmax::T = 10.0
    saveat::T = 0.01
    quadgk_kwargs = (atol=1e-6, order=20)
end

function ηKKa(x, q)
    a = q[1]
    return kink(x + π / 2 + a) + kink(x + π / 2 - a) - 2.0
end

function xspanKKa(q)
    a = q[1]
    return (-π / 2 - abs(a), π / 2 + abs(a))
end

function collective_coordinates(params::KKa)
    @unpack v, tmax, quadgk_kwargs, saveat = params
    q̇₀ = [-v]
    q₀ = [π / 2]

    prob = SecondOrderODEProblem(
        cceq!, q̇₀, q₀, (0.0, tmax), (quadratic, ηKKa, xspanKKa, quadgk_kwargs)
    )
    sol = solve(prob; saveat)
    return Dict("solution" => sol)
end

# Position and Derrick mode

@with_kw struct KKab{T<:Real}
    v::T
    tmax::T = 10.0
    saveat::T = 0.01
    quadgk_kwargs = (atol=1e-6, order=20)
end

function ηKKab(x, q)
    a, b = q
    return kink(b * (x + a) + π / 2) + kink(b * (x - a) + π / 2) - 2.0
end

function xspanKKab(q)
    a, b = q
    return (-π / (2b) - abs(a), π / (2b) + abs(a))
end

function collective_coordinates(params::KKab)
    @unpack v, tmax, quadgk_kwargs, saveat = params
    q̇₀ = [-v, 0.0]
    q₀ = [π / (2 * γ(v)), γ(v)]

    prob = SecondOrderODEProblem(
        cceq!, q̇₀, q₀, (0.0, tmax), (quadratic, ηKKab, xspanKKab, quadgk_kwargs)
    )
    sol = solve(prob; saveat)
    return Dict("solution" => sol)
end
