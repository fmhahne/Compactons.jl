export kink_oscillon_scattering, kink_oscillon_scattering_energies

function getenergies(u, t, integrator)
    φ = @views u[end÷2+1:end]
    ∂ₜφ = @views u[begin:end÷2]

    N = length(φ)
    model, dx = integrator.p

    N₁ = N ÷ 2
    N₂ = N ÷ 2 + round(Int, π / dx)

    E₁ = dx * sum(𝒯(∂ₜφ[i], (φ[i+1] - φ[i-1]) / (2dx)) + model.V(φ[i]) for i ∈ 2:N₁)
    E₂ = dx * sum(𝒯(∂ₜφ[i], (φ[i+1] - φ[i-1]) / (2dx)) + model.V(φ[i]) for i ∈ N₁+1:N₂) - π / 2
    E₃ = dx * sum(𝒯(∂ₜφ[i], (φ[i+1] - φ[i-1]) / (2dx)) + model.V(φ[i]) for i ∈ N₂+1:N-1)

    return [E₁; E₂; E₃]
end

function αₗ(V; v₀=0.0)
    if V ≥ 0
        (-1 + V * (2 + v₀)) / 2
    else
        (1 - V * (2 + v₀)) / 2
    end
end

function α₀(V; v₀=0.0)
    if V ≥ 0
        V
    else
        1 - V
    end
end

function αᵤ(V; v₀=0.0)
    if V ≥ 0
        (1 + V * (2 + v₀)) / 2
    else
        (3 - V * (2 + v₀)) / 2
    end
end

function ab(α, V; v₀=0.0)
    Vc = 1 / (2 + v₀)

    A = (1, 1)
    B = (-1, 0)
    C = (1, 0)
    D = (-1, 1)
    E = (1, -1)
    F = (-1, 2)

    if abs(V) == Vc
        if α ≤ α₀(V; v₀)
            V > 0 ? B : D
        else
            V > 0 ? C : E
        end
    elseif abs(V) < Vc
        if α ≤ α₀(V; v₀)
            V > 0 ? B : D
        elseif α ≤ αᵤ(V; v₀)
            V > 0 ? C : E
        else
            V > 0 ? D : F
        end
    else
        if α ≤ αₗ(V; v₀)
            V > 0 ? A : C
        elseif α ≤ α₀(V; v₀)
            V > 0 ? B : D
        else
            V > 0 ? C : E
        end
    end
end

function x_R(α, V; l=1.0, v₀=0.0)
    a, b = ab(α, V; v₀)
    l * γ(V) / (1 + a * v₀ * V) * ((1 - V^2) * (1 + v₀ * b) + α * (V + v₀ * a))
end

function kink_oscillon_scattering(l, V, α, v₀; dx=1e-3, sampling=10)
    tsave = 0.0:(dx*sampling):10.0

    x = -tsave[end]:dx:tsave[end]
    xsave = x[begin:sampling:end]

    η₀ = kink.(0.0, x) + oscillon.(l * α * γ(V), x .+ x_R(α, V; l, v₀), V; l, v₀)
    ∂ₜη₀ = ∂ₜkink.(0.0, x) + ∂ₜoscillon.(l * α * γ(V), x .+ x_R(α, V; l, v₀), V; l, v₀)

    energies = SavedValues(Float64, Vector{Float64})
    cbenergies = SavingCallback(getenergies, energies; saveat=tsave)

    η, H = producedata(∂ₜη₀, η₀, (quadratic, dx), tsave; callbacks=[cbenergies], dt=0.1dx, sampling)
    E = reduce(hcat, energies.saveval)

    return (x=xsave, t=tsave, η=η, H=H, E₁=E[1, :], E₂=E[2, :], E₃=E[3, :])
end
