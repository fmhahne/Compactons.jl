@with_kw struct KinkAntikink{T<:Real}
    V::T
    dx::T = 1e-3
    dt::T = 0.1dx
    tmax::T = 10.0
    sampling::Int = 10
end

function simulation(parameters::KinkAntikink)
    @unpack V, dx, dt, tmax, sampling = parameters
    tsave = 0.0:(dx * sampling):tmax

    x = (-tmax):dx:tmax
    xsave = x[begin:sampling:end]

    η₀ = kink.(0.0, -abs.(x) .+ π / γ(V), V)
    ∂ₜη₀ = ∂ₜkink.(0, -abs.(x) .+ π / γ(V), V)

    η, H = producedata(quadratic, ∂ₜη₀, η₀, tsave; dx, dt, sampling)
    return Dict("x" => xsave, "t" => tsave, "η" => η, "H" => H)
end
