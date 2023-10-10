using Test
using Compactons

const dx = 1e-3
const x = -2:dx:5
const N = length(x)

@test oscillon.(0.0, x) == oscillon.(1.0, x) == zero(x)
@test oscillon.(-1.0, x, l=2.0) == oscillon.(1.0, x, l=2.0)
@test ∂ₜoscillon.(0.0, x) == ∂ₜoscillon.(1.0, x) <= zero(x)

@test kink.(0.0) == 0.0
@test kink.(π) == 2.0

tsave = 0:1e-2:1.0

ϕ₀ = oscillon.(0.0, x, 0.1)
∂ₜϕ₀ = ∂ₜoscillon.(0.0, x, 0.1)
field1, hamiltonian1 = produce_data(signumgordon, ∂ₜϕ₀, ϕ₀, tsave; dx, sampling=1)

for (i, t) in enumerate(tsave)
    @test field1[:, i] ≈ oscillon.(t, x, 0.1) atol = 1e-3
end

η₀ = kink.(0.0, x, 0.1)
∂ₜη₀ = ∂ₜkink.(0.0, x, 0.1)
field2, hamiltonian2 = produce_data(quadratic, ∂ₜη₀, η₀, tsave; dx, sampling=1)

for (i, t) in enumerate(tsave)
    @test field2[:, i] ≈ kink.(t, x, 0.1) atol = 1e-3
end
