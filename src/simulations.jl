using DifferentialEquations, HDF5

export produce_data, save_data
export kink_antikink_scattering, kink_oscillon_scattering, kink_oscillon_superposition, perturbed_kink

const dx = 8e-4

function produce_data(∂ₜfield₀, field₀, tsave, N, dx, save_idxs; model=quadratic!, hamiltonian=quadratic_hamiltonian)
    saved_hamiltonian = SavedValues(Float64, Vector{Float64})
    callback = SavingCallback(hamiltonian, saved_hamiltonian, saveat=tsave)

    tspan = (tsave[1], tsave[end])
    prob = SecondOrderODEProblem(model, ∂ₜfield₀, field₀, tspan, (N, dx))
    sol = solve(prob, RK4(), adaptive=false, dt=1e-4, saveat=tsave, save_idxs=save_idxs .+ N, callback=callback)

    field = reduce(hcat, sol.u)
    hamiltonian = reduce(hcat, saved_hamiltonian.saveval)
    return field, hamiltonian
end

function save_data(filename, xsave, tsave, field, hamiltonian)
    h5open(filename, "w") do file
        file["x", deflate=3] = collect(xsave)
        file["t", deflate=3] = collect(tsave)
        file["field", deflate=3] = field
        file["hamiltonian", deflate=3] = hamiltonian
    end

    return nothing
end

function xgrid(tsave)
    x = -tsave[end]:dx:tsave[end]
    N = length(x)
    save_idxs = 1:10:N
    xsave = x[save_idxs]

    return x, N, save_idxs, xsave
end

function kink_antikink_scattering(V)
    tsave = ifelse(V < 0.5, 0.0:1e-2:20.0, 0.0:1e-2:10.0)
    x, N, save_idxs, xsave = xgrid(tsave)

    η₀ = kink.(0.0, -abs.(x) .+ π / γ(V), V)
    ∂ₜη₀ = ∂ₜkink.(0, -abs.(x) .+ π / γ(V), V)

    field, hamiltonian = produce_data(∂ₜη₀, η₀, tsave, N, dx, save_idxs)
    return xsave, tsave, field, hamiltonian
end

function kink_oscillon_scattering(l, V, α)
    tsave = ifelse(V < 0.5, 0.0:1e-2:15.0, 0.0:1e-2:10.0)
    x, N, save_idxs, xsave = xgrid(tsave)

    η₀ = kink.(0.0, x .+ π / γ(V), V) + oscillon.(α * l, x, l=l)
    ∂ₜη₀ = ∂ₜkink.(0.0, x .+ π / γ(V), V) + ∂ₜoscillon.(α * l, x, l=l)

    field, hamiltonian = produce_data(∂ₜη₀, η₀, tsave, N, dx, save_idxs)
    return xsave, tsave, field, hamiltonian
end

function kink_oscillon_superposition(l, α)
    tsave = 0.0:1e-2:10.0
    x, N, save_idxs, xsave = xgrid(tsave)

    η₀ = kink.(x .+ π / 2) + oscillon.(α * l, x .+ l / 2, l=l)
    ∂ₜη₀ = ∂ₜoscillon.(α * l, x .+ l / 2, l=l)

    field, hamiltonian = produce_data(∂ₜη₀, η₀, tsave, N, dx, save_idxs)
    return xsave, tsave, field, hamiltonian
end


function perturbed_kink(ϵ)
    tsave = 0.0:1e-2:10.0
    x, N, save_idxs, xsave = xgrid(tsave)

    η₀ = kink.(x / (1.0 + ϵ) .+ π / 2)
    ∂ₜη₀ = zero(x)

    field, hamiltonian = produce_data(∂ₜη₀, η₀, tsave, N, dx, save_idxs)
    return xsave, tsave, field, hamiltonian
end
