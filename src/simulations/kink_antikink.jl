@with_kw struct KinkAntikink{T<:Real}
    v::T
    dx::T = 1e-3
    dt::T = 0.1dx
    tmax::T = 10.0
    sampling::Int = 10
end

function simulation(params::KinkAntikink)
    @unpack v, dx, dt, tmax, sampling = params
    tsave = 0.0:(dx * sampling):tmax

    x = (-tmax):dx:tmax
    xsave = x[begin:sampling:end]

    η₀ = kink.(0.0, -abs.(x) .+ π / γ(v), v)
    ∂ₜη₀ = ∂ₜkink.(0, -abs.(x) .+ π / γ(v), v)

    η, H = produce_data(quadratic, ∂ₜη₀, η₀, tsave; dx, dt, sampling)
    return Dict("x" => xsave, "t" => tsave, "η" => η, "H" => H)
end

@with_kw struct KinkAntikinkMiddle{T<:Real}
    v::T
    dx::T = 1e-3
    dt::T = 0.1dx
    tmax::T = 10.0
    sampling::Int = 10
end

function simulation(params::KinkAntikinkMiddle)
    @unpack v, dx, dt, tmax, sampling = params
    kk_params = KinkAntikink(; v, dx, dt, tmax, sampling)
    @unpack t, η = simulation(kk_params)
    return Dict("t" => t, "η0" => η[end ÷ 2 + 1, :])
end
