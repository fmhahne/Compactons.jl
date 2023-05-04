module KinkOscillonInteraction

using DifferentialEquations
using HypergeometricFunctions: _₂F₁general2 as _₂F₁
using LoopVectorization
using Parameters
using Roots
using SpecialFunctions: gamma as Γ

export γ, boost
export Model
export signumgordon, quadratic, toy, generalizedmodel
export fieldeq!, gethamiltonian, getenergy, 𝒯
export producedata, simulation
export KinkOscillon, NonBPSKink
export kink, ∂ₜkink, ∂ₓkink
export oscillon, ∂ₜoscillon, ∂ₓoscillon
export x_L, x_R, L
export toykink
export generalizedkink, x₀

include("lorentz.jl")
include("models.jl")
include("simulations.jl")
include("solutions/generalizedkink.jl")
include("solutions/kink.jl")
include("solutions/oscillon.jl")
include("solutions/toykink.jl")

end # module
