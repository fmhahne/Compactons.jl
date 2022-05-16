module Compactons

using DifferentialEquations
using LoopVectorization
using Parameters

export γ, boost
export Model
export signumgordon, quadratic, toy
export fieldeq!, gethamiltonian, getenergy, 𝒯
export producedata, simulation
export KinkAntikink, KinkOscillon, NonBPSKink
export kink, ∂ₜkink, ∂ₓkink
export oscillon, ∂ₜoscillon, ∂ₓoscillon
export x_L, x_R, L
export toykink

include("lorentz.jl")
include("models.jl")
include("simulations.jl")
include("solutions/kink.jl")
include("solutions/oscillon.jl")
include("solutions/toykink.jl")

end # module
