module Compactons

using DifferentialEquations
using FiniteDiff
using HypergeometricFunctions: _₂F₁general2 as _₂F₁
using LinearAlgebra
using LoopVectorization
using Parameters
using QuadGK
using Roots
using SpecialFunctions: gamma as Γ

export γ, boost
export Model
export signum_gordon, quadratic, toy, generalized_model
export field_equation!, get_hamiltonian, get_energy, 𝒯
export produce_data, simulation
export KinkAntikink, KinkKink, KinkOscillon, DeformedKink
export kink, ∂ₜkink, ∂ₓkink
export oscillon, ∂ₜoscillon, ∂ₓoscillon
export x_L, x_R, L
export toy_kink
export generalized_kink, x₀
export collective_coordinates
export KinkKinkBorder
export TriangularDelta
export CCDeformedKink
export CCKinkAntikink, ηKAK
export CCKinkKinkNonRel, CCKinkKinkNonRelMode
export CCKinkKinkRel, CCKinkKinkRelMode

include("collective_coordinates.jl")
include("collective_coordinates/deformed_kink.jl")
include("collective_coordinates/kink_antikink.jl")
include("collective_coordinates/kink_kink.jl")
include("lorentz.jl")
include("models.jl")
include("simulations.jl")
include("simulations/deformed_kink.jl")
include("simulations/kink_antikink.jl")
include("simulations/kink_kink.jl")
include("simulations/kink_oscillon.jl")
include("simulations/triangular_delta.jl")
include("solutions/generalized_kink.jl")
include("solutions/kink.jl")
include("solutions/oscillon.jl")
include("solutions/toy_kink.jl")

end # module
