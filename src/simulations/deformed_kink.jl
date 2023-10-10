@with_kw struct DeformedKink{T<:Real}
    ϵ::T
    dx::T = 1e-3
    dt::T = 0.1dx
    tmax::T = 10.0
    sampling::Int = 10
end

function simulation(params::DeformedKink)
    @unpack ϵ, dx, dt, tmax, sampling = params
    tsave = 0.0:(dx * sampling):tmax

    x = (-tmax):dx:tmax
    xsave = x[begin:sampling:end]

    η₀ = kink.(x / (1.0 + ϵ) .+ π / 2)
    ∂ₜη₀ = zero(x)

    η, H = produce_data(quadratic, ∂ₜη₀, η₀, tsave; dx, sampling)
    return Dict("x" => xsave, "t" => tsave, "η" => η, "H" => H)
end
