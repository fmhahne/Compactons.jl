using DifferentialEquations

export producedata
export kink_antikink_scattering, kink_oscillon_scattering, kink_oscillon_superposition, perturbed_kink

const dx = 8e-4

function producedata(∂ₜφ₀, φ₀, p, tsave, save_idxs)
    savedhamiltonian = SavedValues(Float64, Vector{Float64})
    callback = SavingCallback(gethamiltonian, savedhamiltonian, saveat=tsave)

    tspan = (tsave[1], tsave[end])
    prob = SecondOrderODEProblem(fieldeq!, ∂ₜφ₀, φ₀, tspan, p)
    sol = solve(prob, RK4(), adaptive=false, dt=1e-4, saveat=tsave, save_idxs=save_idxs, callback=callback)

    φ = reduce(hcat, sol.u)
    H = reduce(hcat, savedhamiltonian.saveval)
    return φ, H
end

function xgrid(tsave)
    x = -tsave[end]:dx:tsave[end]
    N = length(x)
    save_idxs = N+1:10:2N
    xsave = x[save_idxs.-N]

    return x, N, save_idxs, xsave
end

function kink_antikink_scattering(V)
    tsave = ifelse(V < 0.5, 0.0:1e-2:20.0, 0.0:1e-2:10.0)
    x, N, save_idxs, xsave = xgrid(tsave)

    η₀ = kink.(0.0, -abs.(x) .+ π / γ(V), V)
    ∂ₜη₀ = ∂ₜkink.(0, -abs.(x) .+ π / γ(V), V)

    η, H = producedata(∂ₜη₀, η₀, (quadratic, N, dx), tsave, save_idxs)
    return xsave, tsave, η, H
end

function kink_oscillon_scattering(l, V, α)
    tsave = ifelse(V < 0.5, 0.0:1e-2:15.0, 0.0:1e-2:10.0)
    x, N, save_idxs, xsave = xgrid(tsave)

    η₀ = kink.(0.0, x .+ π / γ(V), V) + oscillon.(α * l, x, l=l)
    ∂ₜη₀ = ∂ₜkink.(0.0, x .+ π / γ(V), V) + ∂ₜoscillon.(α * l, x, l=l)

    η, H = producedata(∂ₜη₀, η₀, (quadratic, N, dx), tsave, save_idxs)
    return xsave, tsave, η, H
end

function kink_oscillon_superposition(l, α)
    tsave = 0.0:1e-2:10.0
    x, N, save_idxs, xsave = xgrid(tsave)

    η₀ = kink.(x .+ π / 2) + oscillon.(α * l, x .+ l / 2, l=l)
    ∂ₜη₀ = ∂ₜoscillon.(α * l, x .+ l / 2, l=l)

    η, H = producedata(∂ₜη₀, η₀, (quadratic, N, dx), tsave, save_idxs)
    return xsave, tsave, η, H
end


function perturbed_kink(ϵ)
    tsave = 0.0:1e-2:10.0
    x, N, save_idxs, xsave = xgrid(tsave)

    η₀ = kink.(x / (1.0 + ϵ) .+ π / 2)
    ∂ₜη₀ = zero(x)

    η, H = producedata(∂ₜη₀, η₀, (quadratic, N, dx), tsave, save_idxs)
    return xsave, tsave, η, H
end
