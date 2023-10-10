module Compactons

using DiffEqPhysics
using DifferentialEquations
using HypergeometricFunctions: _₂F₁general2 as _₂F₁
using LinearAlgebra
using LoopVectorization
using Parameters
using Roots
using SpecialFunctions: gamma as Γ
using TensorOperations

export γ, boost
export Model
export signum_gordon, quadratic, toy, generalized_model
export field_equation!, get_hamiltonian, 𝒯
export produce_data, simulation
export KinkAntikink, KinkKink, KinkOscillon, DeformedKink
export kink, ∂ₜkink, ∂ₓkink
export oscillon, ∂ₜoscillon, ∂ₓoscillon
export x_L, x_R, L
export toykink
export generalized_kink, x₀
export collective_coordinates
export KKa, KKab
export KinkKinkBorder
export TriangularDelta
export moduli_space
export KinkAntikinkNonRelModuliSpace, KinkAntikinkRelModuliSpace
export KinkKinkNonRelModuliSpace
export DeformedKinkModuliSpace

include("collective_coordinates.jl")
include("lorentz.jl")
include("models.jl")
include("moduli_space.jl")
include("moduli_space/deformed_kink.jl")
include("moduli_space/kink_antikink_non_relativistic.jl")
include("moduli_space/kink_antikink_relativistic.jl")
include("moduli_space/kink_kink_non_relativistic.jl")
include("simulations.jl")
include("simulations/deformed_kink.jl")
include("simulations/kink_antikink.jl")
include("simulations/kink_kink.jl")
include("simulations/kink_oscillon.jl")
include("simulations/triangular_delta.jl")
include("solutions/generalized_kink.jl")
include("solutions/kink.jl")
include("solutions/oscillon.jl")
include("solutions/toykink.jl")

end # module
