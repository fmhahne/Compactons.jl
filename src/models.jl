struct Model
    V::Function
    V′::Function
end

signum_gordon = Model(ϕ -> abs(ϕ), ϕ -> sign(ϕ))

quadratic = Model(η -> mod(η, 2) - mod(η, 2)^2 / 2, η -> sign(mod(η, 2)) - mod(η, 2))

toy = Model(η -> abs(mod(η - 1, 2) - 1), η -> sign(mod(η - 1, 2) - sign(mod(η - 1, 2))))

klein_gordon = Model(ϕ -> ϕ^2 / 2, ϕ -> ϕ)

function tanh_gordon(k)
    if k == 0
        return klein_gordon
    end
    return Model(
        ϕ -> exp(-k) * ϕ^2 / 2 + log(cosh(k * ϕ)) / k, ϕ -> exp(-k) * ϕ + tanh(k * ϕ)
    )
end

function generalizedmodel(k)
    return Model(
        η -> (2 - k) / 2 * (1 - abs(mod(η, 2) - 1)^(k + 1)),
        η ->
            -(2 - k) / 2 *
            (k + 1) *
            abs(mod(η, 2) - 1)^k *
            sign(mod(η, 2) - 1) *
            sign(mod(η, 2)),
    )
end

function field_equation!(∂ₜₜϕ, ∂ₜϕ, ϕ, (model, dx), t)
    N = length(ϕ)

    ∂ₜₜϕ[1] = 0
    @tturbo for i in 2:(N - 1)
        ∂ₜₜϕ[i] = (ϕ[i + 1] + ϕ[i - 1] - 2ϕ[i]) / dx^2 - model.V′(ϕ[i])
    end
    ∂ₜₜϕ[N] = 0

    return nothing
end

𝒯(∂ₜϕ, ∂ₓϕ) = (∂ₜϕ^2 + ∂ₓϕ^2) / 2

function get_hamiltonian(u, t, integrator)
    ϕ = @views u[(end ÷ 2 + 1):end]
    ∂ₜϕ = @views u[begin:(end ÷ 2)]

    N = length(ϕ)
    model, dx = integrator.p
    save_idxs = integrator.opts.save_idxs .- N

    H = zero(ϕ)
    for i in intersect(2:(N - 1), save_idxs)
        @inbounds H[i] = 𝒯(∂ₜϕ[i], (ϕ[i + 1] - ϕ[i - 1]) / (2dx)) + model.V(ϕ[i])
    end

    return H[save_idxs]
end
