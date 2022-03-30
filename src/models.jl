using LoopVectorization

export signumgordon!, signumgordon_hamiltonian
export quadratic!, quadratic_hamiltonian
export toy!, toy_hamiltonian

function signumgordon!(∂ₜₜφ, ∂ₜφ, φ, (N, dx), t)
    ∂ₜₜφ[1] = 0.0
    @tturbo for i ∈ 2:N-1
        ∂ₜₜφ[i] = (φ[i+1] + φ[i-1] - 2φ[i]) / dx^2 - sign(φ[i])
    end
    ∂ₜₜφ[N] = 0.0

    nothing
end

function signumgordon_hamiltonian(u, t, integrator)
    N, dx = integrator.p
    save_idxs = integrator.opts.save_idxs .- N

    φ = @views u[N+1:2N]
    ∂ₜφ = @views u[1:N]

    H = zero(φ)
    @inbounds for i ∈ intersect(2:N-1, save_idxs)
        H[i] = ((φ[i+1] - φ[i-1]) / (2dx))^2 / 2 + (∂ₜφ[i])^2 / 2 + abs(φ[i])
    end

    return H[save_idxs]
end

function quadratic!(∂ₜₜη, ∂ₜη, η, (N, dx), t)
    ∂ₜₜη[1] = 0.0
    @tturbo for i ∈ 2:N-1
        ∂ₜₜη[i] = (η[i+1] + η[i-1] - 2η[i]) / dx^2 - sign(mod(η[i], 2)) + mod(η[i], 2)
    end
    ∂ₜₜη[N] = 0.0

    nothing
end

function quadratic_hamiltonian(u, t, integrator)
    N, dx = integrator.p
    save_idxs = integrator.opts.save_idxs .- N

    η = @views u[N+1:2N]
    ∂ₜη = @views u[1:N]

    H = zero(η)
    @inbounds for i ∈ intersect(2:N-1, save_idxs)
        H[i] = ((η[i+1] - η[i-1]) / (2dx))^2 / 2 + (∂ₜη[i])^2 / 2 + abs(mod(η[i], 2)) - mod(η[i], 2)^2 / 2
    end

    return H[save_idxs]
end

function toy!(∂ₜₜη, ∂ₜη, η, (N, dx), t)
    ∂ₜₜη[1] = 0.0
    @tturbo for i ∈ 2:N-1
        ∂ₜₜη[i] = (η[i+1] + η[i-1] - 2η[i]) / dx^2 - sign(mod(η[i] - 1, 2) - 1)
    end
    ∂ₜₜη[N] = 0.0

    nothing
end

function toy_hamiltonian(u, t, integrator)
    N, dx = integrator.p
    save_idxs = integrator.opts.save_idxs .- N

    η = @views u[N+1:2N]
    ∂ₜη = @views u[1:N]

    H = zero(η)
    @inbounds for i ∈ intersect(2:N-1, save_idxs)
        H[i] = ((η[i+1] - η[i-1]) / (2dx))^2 / 2 + (∂ₜη[i])^2 / 2 + abs(mod(η[i] - 1, 2) - 1)
    end

    return H[save_idxs]
end
