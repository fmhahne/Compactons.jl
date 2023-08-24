struct Model
    V::Function
    V′::Function
end

signumgordon = Model(φ -> abs(φ), φ -> sign(φ))

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

function fieldeq!(∂ₜₜφ, ∂ₜφ, φ, (model, dx), t)
    N = length(φ)

    ∂ₜₜφ[1] = 0
    @tturbo for i in 2:(N - 1)
        ∂ₜₜφ[i] = (φ[i + 1] + φ[i - 1] - 2φ[i]) / dx^2 - model.V′(φ[i])
    end
    ∂ₜₜφ[N] = 0

    return nothing
end

𝒯(∂ₜφ, ∂ₓφ) = (∂ₜφ^2 + ∂ₓφ^2) / 2

function gethamiltonian(u, t, integrator)
    φ = @views u[(end ÷ 2 + 1):end]
    ∂ₜφ = @views u[begin:(end ÷ 2)]

    N = length(φ)
    model, dx = integrator.p
    save_idxs = integrator.opts.save_idxs .- N

    H = zero(φ)
    for i in intersect(2:(N - 1), save_idxs)
        @inbounds H[i] = 𝒯(∂ₜφ[i], (φ[i + 1] - φ[i - 1]) / (2dx)) + model.V(φ[i])
    end

    return H[save_idxs]
end

function getenergy(u, t, integrator)
    φ = @views u[(end ÷ 2 + 1):end]
    ∂ₜφ = @views u[begin:(end ÷ 2)]

    N = length(φ)
    model, dx = integrator.p

    return dx *
           sum(𝒯(∂ₜφ[i], (φ[i + 1] - φ[i - 1]) / (2dx)) + model.V(φ[i]) for i in 2:(N - 1))
end
