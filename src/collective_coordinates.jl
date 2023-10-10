function cceq!(q̈, q̇, q, (model, η, xspan, quadgk_kwargs), t)
    ∂η_∂q(x) = FiniteDiff.finite_difference_gradient(q -> η(x, q), q)
    ∂η_∂x(x) = FiniteDiff.finite_difference_derivative(x -> η(x, q), x)
    ∂²η_∂q∂q(x) = FiniteDiff.finite_difference_hessian(q -> η(x, q), q)
    ∂²η_∂q∂x(x) = FiniteDiff.finite_difference_derivative(∂η_∂q, x)

    xi, xf = xspan(q)

    g, _ = quadgk(x -> ∂η_∂q(x) * ∂η_∂q(x)', xi, xf; quadgk_kwargs...)
    g⁻¹ = inv(g)

    f, _ = quadgk(
        x -> ∂²η_∂q∂x(x) * ∂η_∂x(x) + model.V′(η(x, q)) * ∂η_∂q(x), xi, xf; quadgk_kwargs...
    )

    Γ = zeros(eltype(q), length(q), length(q), length(q))
    for (l, j, k) in zip(eachindex(q), eachindex(q), eachindex(q))
        @inbounds Γ[l, j, k], _ = quadgk(
            x -> ∂η_∂q(x)[l] * ∂²η_∂q∂q(x)[j, k], xi, xf; quadgk_kwargs...
        )
    end

    @tensor q̈[i] = -g⁻¹[i, l] * (f[l] + Γ[l, j, k] * q̇[j] * q̇[k])
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
    quadgk_kwargs = (atol=1e-6, order=20)
end

function ηKKa(x, q)
    a = q[1]
    return kink(x + π / 2 + a) + kink(x + π / 2 - a) - 2.0
end

function xspanKKa(q)
    a = q[1]
    return (-π / 2 - a, π / 2 + a)
end

function collective_coordinates(params::KKa)
    @unpack v, tmax, quadgk_kwargs = params
    q̇₀ = [-v]
    q₀ = [π / 2]

    prob = SecondOrderODEProblem(
        cceq!, q̇₀, q₀, (0.0, tmax), (quadratic, ηKKa, xspanKKa, quadgk_kwargs)
    )
    return Dict("solution" => solve(prob; saveat=1e-2))
end

# Position and Derrick mode

@with_kw struct KKab{T<:Real}
    v::T
    tmax::T = 10.0
    quadgk_kwargs = (atol=1e-6, order=20)
end

function ηKKab(x, q)
    a, b = q
    return kink(b * (x + a) + π / 2) + kink(b * (x - a) + π / 2) - 2.0
end

function xspanKKab(q)
    a, b = q
    return (-π / (2b) - a, π / (2b) + a)
end

function collective_coordinates(params::KKab)
    @unpack v, tmax, quadgk_kwargs = params
    q̇₀ = [-v, 0.0]
    q₀ = [π / (2 * γ(v)), γ(v)]

    prob = SecondOrderODEProblem(
        cceq!, q̇₀, q₀, (0.0, tmax), (quadratic, ηKKab, xspanKKab, quadgk_kwargs)
    )
    return Dict("solution" => solve(prob))
end
