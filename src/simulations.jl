using DifferentialEquations

export producedata
export kink_antikink_scattering, kink_oscillon_scattering, kink_oscillon_superposition, perturbed_kink

function producedata(∂ₜφ₀, φ₀, (model, N, dx), tsave; sampling=10, callbacks=[])
    savedhamiltonian = SavedValues(Float64, Vector{Float64})
    cbhamiltonian = SavingCallback(gethamiltonian, savedhamiltonian, saveat=tsave)

    callback = CallbackSet(cbhamiltonian, callbacks...)

    tspan = (tsave[begin], tsave[end])
    prob = SecondOrderODEProblem(fieldeq!, ∂ₜφ₀, φ₀, tspan, (model, N, dx))
    sol = solve(prob, RK4(); adaptive=false, dt=0.125dx, saveat=tsave, save_idxs=N+1:sampling:2N, callback=callback)

    φ = reduce(hcat, sol.u)
    H = reduce(hcat, savedhamiltonian.saveval)
    return φ, H
end

function kink_antikink_scattering(V; dx=8e-4, sampling=10)
    tsave = ifelse(V < 0.5, 0.0:1e-2:20.0, 0.0:1e-2:10.0)

    x = -tsave[end]:dx:tsave[end]
    xsave = x[begin:sampling:end]
    N = length(x)

    η₀ = kink.(0.0, -abs.(x) .+ π / γ(V), V)
    ∂ₜη₀ = ∂ₜkink.(0, -abs.(x) .+ π / γ(V), V)

    η, H = producedata(∂ₜη₀, η₀, (quadratic, N, dx), tsave; sampling=sampling)
    return xsave, tsave, η, H
end

function kink_oscillon_scattering(l, V, α; dx=8e-4, sampling=10)
    tsave = ifelse(V < 0.5, 0.0:1e-2:15.0, 0.0:1e-2:10.0)

    x = -tsave[end]:dx:tsave[end]
    xsave = x[begin:sampling:end]
    N = length(x)

    η₀ = kink.(0.0, x .+ π / γ(V), V) + oscillon.(α * l, x, l=l)
    ∂ₜη₀ = ∂ₜkink.(0.0, x .+ π / γ(V), V) + ∂ₜoscillon.(α * l, x, l=l)

    η, H = producedata(∂ₜη₀, η₀, (quadratic, N, dx), tsave; sampling=sampling)
    return xsave, tsave, η, H
end

function kink_oscillon_superposition(l, α; dx=8e-4, sampling=10)
    tsave = 0.0:1e-2:10.0

    x = -tsave[end]:dx:tsave[end]
    xsave = x[begin:sampling:end]
    N = length(x)

    η₀ = kink.(x .+ π / 2) + oscillon.(α * l, x .+ l / 2, l=l)
    ∂ₜη₀ = ∂ₜoscillon.(α * l, x .+ l / 2, l=l)

    η, H = producedata(∂ₜη₀, η₀, (quadratic, N, dx), tsave; sampling=sampling)
    return xsave, tsave, η, H
end


function perturbed_kink(ϵ; dx=8e-4, sampling=10)
    tsave = 0.0:1e-2:10.0

    x = -tsave[end]:dx:tsave[end]
    xsave = x[begin:sampling:end]
    N = length(x)

    η₀ = kink.(x / (1.0 + ϵ) .+ π / 2)
    ∂ₜη₀ = zero(x)

    η, H = producedata(∂ₜη₀, η₀, (quadratic, N, dx), tsave; sampling=sampling)
    return xsave, tsave, η, H
end
