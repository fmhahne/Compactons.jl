using DifferentialEquations

export producedata
export kink_antikink_scattering, kink_oscillon_scattering, kink_oscillon_superposition, perturbed_kink

function producedata(model, ∂ₜφ₀, φ₀, tsave; dx, dt=dx / 10, sampling=10, callbacks=[])
    savedhamiltonian = SavedValues(Float64, Vector{Float64})
    cbhamiltonian = SavingCallback(gethamiltonian, savedhamiltonian, saveat=tsave)

    callback = CallbackSet(cbhamiltonian, callbacks...)

    N = length(φ₀)
    tspan = (tsave[begin], tsave[end])

    prob = SecondOrderODEProblem(fieldeq!, ∂ₜφ₀, φ₀, tspan, (model, dx))
    sol = solve(prob, RK4(); adaptive=false, dt=dt, saveat=tsave, save_idxs=N+1:sampling:2N, callback=callback)

    φ = reduce(hcat, sol.u)
    H = reduce(hcat, savedhamiltonian.saveval)
    return φ, H
end

function kink_antikink_scattering(V; dx=8e-4, sampling=10)
    tsave = ifelse(V < 0.5, 0.0:(dx*sampling):20.0, 0.0:(dx*sampling):10.0)

    x = -tsave[end]:dx:tsave[end]
    xsave = x[begin:sampling:end]

    η₀ = kink.(0.0, -abs.(x) .+ π / γ(V), V)
    ∂ₜη₀ = ∂ₜkink.(0, -abs.(x) .+ π / γ(V), V)

    η, H = producedata(quadratic, ∂ₜη₀, η₀, tsave; dx, sampling=sampling, dt=1e-4)
    return xsave, tsave, η, H
end

function kink_oscillon_superposition(l, α; dx=1e-3, sampling=10)
    tsave = 0.0:(dx*sampling):10.0

    x = -tsave[end]:dx:tsave[end]
    xsave = x[begin:sampling:end]

    η₀ = kink.(x .+ π / 2) + oscillon.(α * l, x .+ l / 2, l=l)
    ∂ₜη₀ = ∂ₜoscillon.(α * l, x .+ l / 2, l=l)

    η, H = producedata(quadratic, ∂ₜη₀, η₀, tsave; dx, sampling)
    return xsave, tsave, η, H
end

function perturbed_kink(ϵ; dx=1e-3, sampling=10)
    tsave = 0.0:(dx*sampling):10.0

    x = -tsave[end]:dx:tsave[end]
    xsave = x[begin:sampling:end]

    η₀ = kink.(x / (1.0 + ϵ) .+ π / 2)
    ∂ₜη₀ = zero(x)

    η, H = producedata(quadratic, ∂ₜη₀, η₀, tsave; dx, sampling)
    return xsave, tsave, η, H
end

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

function kink_oscillon_scattering(l, V, α, v₀; dx=1e-3, sampling=10)
    tsave = 0.0:(dx*sampling):10.0

    x = -tsave[end]:dx:tsave[end]
    xsave = x[begin:sampling:end]

    η₀ = kink.(0.0, x) + oscillon.(l * α * γ(V), x .+ x_R(α, V; l, v₀), V; l, v₀)
    ∂ₜη₀ = ∂ₜkink.(0.0, x) + ∂ₜoscillon.(l * α * γ(V), x .+ x_R(α, V; l, v₀), V; l, v₀)

    energies = SavedValues(Float64, Vector{Float64})
    cbenergies = SavingCallback(getenergies, energies; saveat=tsave)

    η, H = producedata(quadratic, ∂ₜη₀, η₀, tsave; dx, sampling)
    E = reduce(hcat, energies.saveval)

    return (x=xsave, t=tsave, η=η, H=H, E₁=E[1, :], E₂=E[2, :], E₃=E[3, :])
end
