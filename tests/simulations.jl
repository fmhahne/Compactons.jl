using Test
using Compactons

dx = 1e-3
x = -2:dx:5
N = length(x)

@test oscillon.(0.0, x) == oscillon.(1.0, x) == zero(x)
@test oscillon.(-1.0, x, l=2.0) == oscillon.(1.0, x, l=2.0)
@test ∂ₜoscillon.(0.0, x) == ∂ₜoscillon.(1.0, x) <= zero(x)

@test kink.(0.0) == 0.0
@test kink.(π) == 2.0

tsave = 0:1e-2:1.0

φ₀ = oscillon.(0.0, x, 0.1)
∂ₜφ₀ = ∂ₜoscillon.(0.0, x, 0.1)
field1, hamiltonian1 = produce_data(∂ₜφ₀, φ₀, tsave, N, dx, 1:N; model=signumgordon!, hamiltonian=signumgordon_hamiltonian)

for (i, t) ∈ enumerate(tsave)
    @test field1[:, i] ≈ oscillon.(t, x, 0.1) atol = 1e-3
end

η₀ = kink.(0.0, x, 0.1)
∂ₜη₀ = ∂ₜkink.(0.0, x, 0.1)
field2, hamiltonian2 = produce_data(∂ₜη₀, η₀, tsave, N, dx, 1:N; model=quadratic!, hamiltonian=quadratic_hamiltonian)

for (i, t) ∈ enumerate(tsave)
    @test field2[:, i] ≈ kink.(t, x, 0.1) atol = 1e-3
end
