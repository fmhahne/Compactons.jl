module Compactons

using DifferentialEquations
using HypergeometricFunctions: _₂F₁general2 as _₂F₁
using LoopVectorization
using Parameters
using Roots
using SpecialFunctions: gamma as Γ
using TensorOperations

export γ, boost
export Model
export signumgordon, quadratic, toy, generalizedmodel
export fieldeq!, gethamiltonian, getenergy, 𝒯
export producedata, simulation
export KinkAntikink, KinkKink, KinkOscillon, DeformedKink
export kink, ∂ₜkink, ∂ₓkink
export oscillon, ∂ₜoscillon, ∂ₓoscillon
export x_L, x_R, L
export toykink
export generalizedkink, x₀
export collectivecoordinates
export KKa, KKab
export KinkKinkBorder
export TriangularDelta

include("lorentz.jl")
include("models.jl")
include("simulations.jl")
include("simulations/deformed_kink.jl")
include("simulations/kink_antikink.jl")
include("simulations/kink_kink.jl")
include("simulations/kink_oscillon.jl")
include("simulations/triangular_delta.jl")
include("collectivecoordinates.jl")
include("solutions/generalizedkink.jl")
include("solutions/kink.jl")
include("solutions/oscillon.jl")
include("solutions/toykink.jl")

end # module
