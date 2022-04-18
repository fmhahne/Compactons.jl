using LoopVectorization

export Model
export fieldeq!, hamiltonian
export signumgordon, quadratic, toy

struct Model
    V::Function
    V′::Function
end

signumgordon = Model(
    φ -> abs(φ),
    φ -> sign(φ)
)

quadratic = Model(
    η -> abs(mod(η, 2)) - mod(η, 2)^2 / 2,
    η -> sign(mod(η, 2)) - mod(η, 2)
)

toy = Model(
    η -> abs(mod(η - 1, 2) - 1),
    η -> sign(mod(η - 1, 2) - sign(mod(η - 1, 2)))
)

function fieldeq!(∂ₜₜφ, ∂ₜφ, φ, (model, N, dx), t)
    ∂ₜₜφ[1] = 0
    @tturbo for i ∈ 2:N-1
        ∂ₜₜφ[i] = (φ[i+1] + φ[i-1] - 2φ[i]) / dx^2 - model.V′(φ[i])
    end
    ∂ₜₜφ[N] = 0

    nothing
end

function hamiltonian(u, t, integrator)
    model, N, dx = integrator.p
    save_idxs = integrator.opts.save_idxs .- N

    φ = @views u[N+1:2N]
    ∂ₜφ = @views u[1:N]

    H = zero(φ)
    @inbounds for i ∈ intersect(2:N-1, save_idxs)
        H[i] = ((φ[i+1] - φ[i-1]) / (2dx))^2 / 2 + (∂ₜφ[i])^2 / 2 + model.V(φ[i])
    end

    return H[save_idxs]
end
