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
export KinkAntikink, KinkAntikinkMiddle, KinkAntikinkRadiationEnergy
export KinkKink, KinkKinkBorder
export KinkOscillon
export DeformedKink
export SymOscillonScattering
export TriangularDelta
export kink, ∂ₜkink, ∂ₓkink
export oscillon, ∂ₜoscillon, ∂ₓoscillon
export x_L, x_R, L
export toy_kink
export generalized_kink, x₀
export collective_coordinates, potential, metric
export CCDeformedKink
export CCKinkAntikink, ηKAK
export CCKinkKinkNonRel, xspanKK_non_rel, ηKK_non_rel, ηKK_non_rel_mode
export CCKinkKinkNonRelMode, xspanKK_rel, ηKK_rel, ηKK_rel_mode
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
include("simulations/oscillon_scattering.jl")
include("simulations/triangular_delta.jl")
include("solutions/generalized_kink.jl")
include("solutions/kink.jl")
include("solutions/oscillon.jl")
include("solutions/toy_kink.jl")

end # module
