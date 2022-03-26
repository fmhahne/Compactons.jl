using Printf
using DifferentialEquations, HDF5
using Compactons

mkpath("data")

const dx = 8e-4
const x = -12.5:dx:12.5
const N = length(x)

Threads.@threads for V ∈ 0.00:1e-2:0.99
    id = @sprintf "V=%.2f" V
    println("Starting simulation for " * id)

    η₀ = kink.(0.0, -abs.(x) .+ π / γ(V), V)
    ∂ₜη₀ = ∂ₜkink.(0, -abs.(x) .+ π / γ(V), V)

    if V < 0.5
        tsave = 0.0:1e-2:20.0
    else
        tsave = 0.0:1e-2:10.0
    end
    tspan = (tsave[1], tsave[end])

    hamiltonian = SavedValues(Float64, Vector{Float64})
    callback = SavingCallback(quadraticHamiltonian, hamiltonian, saveat=tsave)

    prob = SecondOrderODEProblem(quadratic!, ∂ₜη₀, η₀, tspan, (N, dx))
    sol = solve(prob, RK4(), adaptive=false, dt=1e-4, saveat=tsave, save_idxs=N+1:10:2N, callback=callback)

    if sol.retcode == :Success
        h5open("data/" * id * ".h5", "w") do f
            f["field", deflate=1] = reduce(hcat, sol.u)
            f["hamiltonian", deflate=1] = reduce(hcat, hamiltonian.saveval)
            f["t", deflate=1] = collect(tsave)
            f["x", deflate=1] = collect(x[1:10:end])
        end
    end
end
