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

@with_kw struct KinkAntikinkRadiationEnergy{T<:Real}
    v::T
    dx::T = 1e-3
    dt::T = 0.1dx
    tmax::T = 10.0
    sampling::Int = 10
end

function get_radiation_energy(u, t, integrator)
    η = @views u[(end ÷ 2 + 1):end]
    ∂ₜη = @views u[begin:(end ÷ 2)]

    N = length(η)
    model, dx = integrator.p

    N0 = N ÷ 2 + round(Int, π / dx)

    E_rad =
        dx * sum(
            𝒯(∂ₜη[i], (η[i + 1] - η[i - 1]) / (2 * dx)) + model.V(η[i]) for
            i in (N0 + 1):(N - 1)
        )

    return E_rad
end

function simulation(params::KinkAntikinkRadiationEnergy)
    @unpack v, dx, dt, tmax, sampling = params
    tsave = 0.0:(dx * sampling):tmax

    x = (-tmax):dx:tmax

    η₀ = kink.(0.0, -abs.(x) .+ π / γ(v), v)
    ∂ₜη₀ = ∂ₜkink.(0, -abs.(x) .+ π / γ(v), v)

    E_rad = SavedValues(Float64, Float64)
    cb_E_rad = SavingCallback(get_radiation_energy, E_rad; saveat=tsave)

    produce_data(quadratic, ∂ₜη₀, η₀, tsave; dx, dt, sampling, callbacks=[cb_E_rad])
    return Dict("t" => tsave, "E_rad" => E_rad.saveval)
end
