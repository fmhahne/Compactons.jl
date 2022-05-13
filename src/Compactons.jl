module Compactons

using DifferentialEquations
using Parameters

export γ, boost
export producedata, simulation
export KinkAntikink, KinkOscillon, NonBPSKink
export kink, ∂ₜkink, ∂ₓkink
export oscillon, ∂ₜoscillon, ∂ₓoscillon
export x_L, x_R, L

include("lorentz.jl")
include("models.jl")
include("simulations.jl")
include("solutions/kink.jl")
include("solutions/oscillon.jl")

end # module
