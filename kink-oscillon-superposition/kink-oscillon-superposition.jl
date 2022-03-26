using Printf
using DifferentialEquations, HDF5
using Compactons

mkpath("data")

const dx = 8e-4

for l ∈ 0.5:0.5:3.0
    for α ∈ 0.00:0.25:0.75
        id = @sprintf "l=%.2f,alpha=%.2f" l α
        println("Starting simulation for " * id)

        tsave = 0.0:1e-2:10.0
        tspan = (tsave[1], tsave[end])

        x = -tsave[end]:dx:tsave[end]
        N = length(x)

        η₀ = kink.(x .+ π / 2) + oscillon.(α * l, x .+ l / 2, l=l)
        ∂ₜη₀ = ∂ₜoscillon.(α * l, x .+ l / 2, l=l)

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
end
