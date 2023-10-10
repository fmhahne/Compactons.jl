@with_kw struct KinkOscillon{T<:Real}
    l::T
    V::T
    α::T
    v₀::T
    x₀::T = x_R(α, V; l, v₀)
    model::Model = quadratic
end

function getenergies(u, t, integrator)
    ϕ = @views u[(end ÷ 2 + 1):end]
    ∂ₜϕ = @views u[begin:(end ÷ 2)]

    N = length(ϕ)
    model, dx = integrator.p

    N₁ = N ÷ 2
    N₂ = N ÷ 2 + round(Int, π / dx)

    E₁ = dx * sum(𝒯(∂ₜϕ[i], (ϕ[i + 1] - ϕ[i - 1]) / (2dx)) + model.V(ϕ[i]) for i in 2:N₁)
    E₂ =
        dx *
        sum(𝒯(∂ₜϕ[i], (ϕ[i + 1] - ϕ[i - 1]) / (2dx)) + model.V(ϕ[i]) for i in (N₁ + 1):N₂) -
        π / 2
    E₃ =
        dx * sum(
            𝒯(∂ₜϕ[i], (ϕ[i + 1] - ϕ[i - 1]) / (2dx)) + model.V(ϕ[i]) for
            i in (N₂ + 1):(N - 1)
        )

    return [E₁; E₂; E₃]
end

function simulation(parameters::KinkOscillon; dx=1e-3, sampling=10)
    @unpack l, V, α, v₀, x₀, model = parameters
    tsave = 0.0:(dx * sampling):10.0

    x = (-tsave[end]):dx:tsave[end]
    xsave = x[begin:sampling:end]

    η₀ = oscillon.(l * α * γ(V), x .+ x₀, V; l, v₀)
    ∂ₜη₀ = ∂ₜoscillon.(l * α * γ(V), x .+ x₀, V; l, v₀)

    η₀ += if model == quadratic
        kink.(0, x)
    elseif model == toy
        toykink.(0, x)
    else
        error("Kink not implemented")
    end

    energies = SavedValues(Float64, Vector{Float64})
    cbenergies = SavingCallback(getenergies, energies; saveat=tsave)

    η, H = produce_data(model, ∂ₜη₀, η₀, tsave; dx, sampling, callbacks=[cbenergies])
    E = reduce(hcat, energies.saveval)

    return Dict(
        "x" => xsave,
        "t" => tsave,
        "η" => η,
        "H" => H,
        "E₁" => E[1, :],
        "E₂" => E[2, :],
        "E₃" => E[3, :],
    )
end
