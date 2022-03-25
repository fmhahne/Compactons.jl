using Printf
using DifferentialEquations, HDF5

include("../AnalyticalSolutions.jl")
include("../VShapeModels.jl")

mkpath("data")

const dx = 8e-4
const x = -10.0:dx:10.0
const N = length(x)

const tsave = 0.0:1e-2:10.0
const tspan = (tsave[1], tsave[end])

for ϵ ∈ -0.50:0.01:0.50
    id = @sprintf "eps=%.2f" ϵ
    println("Starting simulation for " * id)

    η₀ = kink.(x / (1.0 + ϵ) .+ π / 2)
    ∂ₜη₀ = zero(x)

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
