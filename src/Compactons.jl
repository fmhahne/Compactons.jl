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
export moduli_space, KinkAntikinkNonRelModuliSpace, KinkAntikinkRelModuliSpace
export DeformedKinkModuliSpace

include("collectivecoordinates.jl")
include("lorentz.jl")
include("models.jl")
include("moduli_space.jl")
include("moduli_space/deformed_kink.jl")
include("moduli_space/kink_antikink_non_relativistic.jl")
include("moduli_space/kink_antikink_relativistic.jl")
include("simulations.jl")
include("simulations/deformed_kink.jl")
include("simulations/kink_antikink.jl")
include("simulations/kink_kink.jl")
include("simulations/kink_oscillon.jl")
include("simulations/triangular_delta.jl")
include("solutions/generalizedkink.jl")
include("solutions/kink.jl")
include("solutions/oscillon.jl")
include("solutions/toykink.jl")

end # module
